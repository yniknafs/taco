'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
from operator import itemgetter
from multiprocessing import Process, JoinableQueue, Value, Lock

from gtf import GTF
from base import Strand
from transfrag import Transfrag
from locus import Locus
from cpathfinder import find_paths
from path_graph import PathGraphFactory, reconstruct_path

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class Config(object):
    def __init__(self, **kwargs):
        self.unresolved_bg_fhs = []
        self.resolved_bg_fhs = []
        self.splice_bed_fh = None
        self.expr_h5fh = None
        self.splice_graph_gtf_fh = None
        self.path_graph_stats_fh = None
        self.assembly_loss_gtf_fh = None
        self.assembly_gtf_fh = None
        self.assembly_bed_fh = None
        # shared memory values (for parallelism)
        self.locus_id_value_obj = LockValue(1)
        self.gene_id_value_obj = LockValue(1)
        self.tss_id_value_obj = LockValue(1)
        self.t_id_value_obj = LockValue(1)
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    @staticmethod
    def defaults():
        self = Config()
        self.guided_strand = False
        self.guided_ends = False
        self.guided_assembly = False
        self.change_point = True
        self.change_point_pvalue = 0.05
        self.change_point_fold_change = 0.8
        self.change_point_trim = True
        self.path_graph_kmax = 0
        self.path_graph_loss_threshold = 0.10
        self.path_frac = 0
        self.max_paths = 0
        self.isoform_frac = 0
        self.max_isoforms = 0
        return self


class Isoform(object):
    __slots__ = ('expr', 'path', 'rel_frac', 'abs_frac',
                 'gene_id', 'tss_id')

    def __init__(self, path=None, expr=0.0, rel_frac=1.0, abs_frac=1.0):
        self.path = path
        self.expr = expr
        self.rel_frac = rel_frac
        self.abs_frac = abs_frac
        self.gene_id = -1
        self.tss_id = -1


class Cluster(object):

    def __init__(self):
        self._id = 0
        self.expr = 1e-10
        self.nodes = set()
        self.paths = []

    def __len__(self):
        return len(self.paths)

    def _add(self, path, expr):
        self.expr += expr
        self.nodes.update(path)
        self.paths.append((path, expr))

    def iterpaths(self):
        if len(self.paths) == 0:
            return
        best_expr = self.paths[0][1]
        for path, expr in self.paths:
            rel_frac = 0.0 if best_expr == 0.0 else expr / best_expr
            abs_frac = 0.0 if self.expr == 0.0 else expr / self.expr
            yield path, expr, rel_frac, abs_frac

    def merge(self, *others):
        for other in others:
            self.expr += other.expr
            self.nodes.update(other.nodes)
            self.paths.extend(other.paths)
        self.paths.sort(key=itemgetter(1), reverse=True)

    @staticmethod
    def build(paths, min_frac=0.0):
        filtered = []
        clusters = {}
        _id = 0
        for i in xrange(len(paths)):
            path, expr = paths[i]
            # does path overlap existing clusters?
            matches = [c for c in clusters.itervalues()
                       if not c.nodes.isdisjoint(path)]
            if len(matches) == 0:
                # make new cluster
                c = Cluster()
                c._id = _id
                _id += 1
                clusters[c._id] = c
            else:
                # check frac in all clusters
                discard = False
                for c in clusters.itervalues():
                    best_expr = c.paths[0][1]
                    rel_frac = 0.0 if best_expr == 0.0 else expr / best_expr
                    if rel_frac < min_frac:
                        discard = True
                        break
                if discard:
                    filtered.append((path, expr))
                    continue
                # merge clusters
                c = matches[0]
                c.merge(*matches[1:])
                for c2 in matches[1:]:
                    del clusters[c2._id]
            # add to cluster
            c._add(path, expr)
        return clusters.values(), filtered


class LockValue(object):
    def __init__(self, initval=0):
        self.val = Value('L', initval)
        self.lock = Lock()

    def next(self):
        with self.lock:
            cur_val = self.val.value
            self.val.value += 1
            return cur_val


def get_gtf_features(chrom, strand, exons, locus_id, gene_id, tss_id,
                     transcript_id, expr, rel_frac, abs_frac):
    tx_start = exons[0].start
    tx_end = exons[-1].end
    strand_str = Strand.to_gtf(strand)
    attr_dict = {'locus_id': locus_id,
                 'gene_id': gene_id,
                 'tss_id': tss_id,
                 'transcript_id': transcript_id}
    f = GTF.Feature()
    f.seqid = chrom
    f.source = 'taco'
    f.feature = 'transcript'
    f.start = tx_start
    f.end = tx_end
    f.score = int(round(1000.0 * rel_frac))
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'expr': '%.3f' % expr,
               'rel_frac': '%.5f' % rel_frac,
               'abs_frac': '%.5f' % abs_frac}
    f.attrs.update(attr_dict)
    yield f
    for e in exons:
        f = GTF.Feature()
        f.seqid = chrom
        f.source = 'taco'
        f.feature = 'exon'
        f.start = e.start
        f.end = e.end
        f.score = int(round(1000.0 * rel_frac))
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {}
        f.attrs.update(attr_dict)
        yield f


def write_bed(chrom, name, strand, score, exons):
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)
    # make bed fields
    fields = [chrom,
              str(tx_start),
              str(tx_end),
              str(name),
              str(score),
              Strand.to_gtf(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str, block_sizes)) + ',',
              ','.join(map(str, block_starts)) + ',']
    return fields


def assign_ids(isoforms, strand, gene_id_value_obj, tss_id_value_obj):
    # map tss positions to unique ids
    tss_pos_id_map = {}
    gene_id = gene_id_value_obj.next()
    for isoform in isoforms:
        start = isoform.path[0].start
        end = isoform.path[-1].end
        # map TSS positions to IDs
        tss_pos = end if strand == Strand.NEG else start
        if tss_pos not in tss_pos_id_map:
            tss_id = tss_id_value_obj.next()
            tss_pos_id_map[tss_pos] = tss_id
        else:
            tss_id = tss_pos_id_map[tss_pos]
        isoform.gene_id = gene_id
        isoform.tss_id = tss_id


def assemble_isoforms(sgraph, config):
    # read in transfrag paths
    pgf = PathGraphFactory(sgraph)
    K, k = pgf.create_optimal(kmax=config.path_graph_kmax,
                              loss_threshold=config.path_graph_loss_threshold,
                              stats_fh=config.path_graph_stats_fh)
    if K is None or len(K) == 0:
        return []
    # smooth kmer graph
    K.apply_smoothing()

    logging.debug('%s:%d-%d[%s] finding paths in k=%d graph '
                  '(%d kmers) source_expr=%f' %
                  (sgraph.chrom, sgraph.start, sgraph.end,
                   Strand.to_gtf(sgraph.strand), k, len(K),
                   K.exprs[K.SOURCE_ID]))
    paths = []
    for kmer_path, expr in find_paths(K, config.path_frac, config.max_paths):
        path = reconstruct_path(kmer_path, K, sgraph)
        logging.debug("\texpr=%f length=%d" % (expr, len(path)))
        paths.append((path, expr))
    # build gene clusters
    clusters, filtered = Cluster.build(paths, min_frac=config.isoform_frac)
    logging.debug('\tclusters: %d filtered: %d' %
                  (len(clusters), len(filtered)))
    gene_isoforms = []
    for cluster in clusters:
        isoforms = []
        for path, expr, rel_frac, abs_frac in cluster.iterpaths():
            isoforms.append(Isoform(path=path, expr=expr, rel_frac=rel_frac,
                                    abs_frac=abs_frac))
        # apply max isoforms limit (per cluster)
        if config.max_isoforms > 0:
            isoforms = isoforms[:config.max_isoforms]
        gene_isoforms.append(isoforms)
    return gene_isoforms


def assemble_gene(sgraph, locus_id_str, config):
    logging.debug('%s:%d-%d[%s] nodes=%d' %
                  (sgraph.chrom, sgraph.start, sgraph.end,
                   Strand.to_gtf(sgraph.strand), len(sgraph.G)))
    # output splice graph node data
    for f in sgraph.get_node_gtf():
        print >>config.splice_graph_gtf_fh, str(f)

    if config.change_point:
        # detect change points
        changepts = sgraph.detect_change_points(
            pval=config.change_point_pvalue,
            fc_cutoff=config.change_point_fold_change)
        logging.debug('%s:%d-%d[%s] change points: %d' %
                      (sgraph.chrom, sgraph.start, sgraph.end,
                       Strand.to_gtf(sgraph.strand), len(changepts)))
        for cp in changepts:
            sgraph.apply_change_point(cp, config.change_point_trim)
            # output splice graph change points
            for f in sgraph.get_change_point_gtf(cp):
                print >>config.splice_graph_gtf_fh, str(f)
        # must recreate splice graph after finding change points
        if len(changepts) > 0:
            sgraph.recreate()

    # run isoform path finding algorithm, filter and group into genes
    for gene_isoforms in assemble_isoforms(sgraph, config):
        # assign gene_id and tss_id
        assign_ids(gene_isoforms, sgraph.strand, config.gene_id_value_obj,
                   config.tss_id_value_obj)
        # write output
        for isoform in gene_isoforms:
            # assign transcript id
            t_id = config.t_id_value_obj.next()
            # get strings for each id
            t_id_str = "TU%d" % t_id
            tss_id_str = "TSS%d" % (isoform.tss_id)
            gene_id_str = "G%d" % (isoform.gene_id)
            # write to GTF
            for f in get_gtf_features(chrom=sgraph.chrom,
                                      strand=sgraph.strand,
                                      exons=isoform.path,
                                      locus_id=locus_id_str,
                                      gene_id=gene_id_str,
                                      tss_id=tss_id_str,
                                      transcript_id=t_id_str,
                                      expr=isoform.expr,
                                      rel_frac=isoform.rel_frac,
                                      abs_frac=isoform.abs_frac):
                print >>config.assembly_gtf_fh, str(f)
            # write to BED
            name = "%s|%s(%.1f)" % (gene_id_str, t_id_str, isoform.expr)
            fields = write_bed(sgraph.chrom, name, sgraph.strand,
                               int(round(1000.0 * isoform.rel_frac)),
                               isoform.path)
            print >>config.assembly_bed_fh, '\t'.join(fields)


def assemble_locus(gtf_lines, config):
    t_dict = Transfrag.parse_gtf(gtf_lines)
    locus = Locus.create(t_dict.values(),
                         config.guided_strand,
                         config.guided_ends,
                         config.guided_assembly)
    logging.debug('\t%d transfrags (+: %d, -: %d, .: %d)' %
                  (len(t_dict),
                   len(locus.get_transfrags(Strand.POS)),
                   len(locus.get_transfrags(Strand.NEG)),
                   len(locus.get_transfrags(Strand.NA))))
    # write raw bedgraph files
    locus.write_bedgraph(config.unresolved_bg_fhs)
    # resolve unstranded transcripts
    num_resolved = locus.impute_unknown_strands()
    if num_resolved > 0:
        logging.debug('\t%d resolved (+: %d, -: %d, .: %d)' %
                      (num_resolved,
                       len(locus.get_transfrags(Strand.POS)),
                       len(locus.get_transfrags(Strand.NEG)),
                       len(locus.get_transfrags(Strand.NA))))
    # write bedgraph files after strand resolved
    locus.write_bedgraph(config.resolved_bg_fhs)
    # write splice junctions
    locus.write_splice_bed(config.splice_bed_fh)
    # convert to stranded locus objects
    locus_id_str = "L%d" % (config.locus_id_value_obj.next())
    for sgraph in locus.create_splice_graphs():
        assemble_gene(sgraph, locus_id_str, config)


def assemble(**kwargs):
    '''
    kwargs: dict containing arguments and input/output file locations

    Configuration attributes:
    - guided_strand
    - guided_ends
    - guided_assembly
    - change_point
    - change_point_pvalue
    - change_point_fold_change
    - change_point_trim
    - path_graph_kmax
    - path_graph_loss_threshold
    - path_frac
    - max_paths
    - isoform_frac
    - max_isoforms

    Input file attributes:
    - transfrags_gtf_file
    - chrom_sizes_file

    Output file attributes:
    - unresolved_bg_files
    - resolved_bg_files
    - splice_bed_file
    - expr_h5_file
    - splice_graph_gtf_file
    - path_graph_stats_file
    - assembly_loss_gtf_file
    - assembly_gtf_file
    - assembly_bed_file
    '''
    config = Config(**kwargs)
    # setup bedgraph output files
    for s, filename in config.unresolved_bg_files:
        config.unresolved_bg_fhs.append(open(filename, 'w'))
    for s, filename in config.resolved_bg_files:
        config.resolved_bg_fhs.append(open(filename, 'w'))
    # setup junction bed file
    config.splice_bed_fh = Locus.open_splice_bed(config.splice_bed_file)
    # splice graph gtf file
    config.splice_graph_gtf_fh = open(config.splice_graph_gtf_file, 'w')
    # path graph stats file
    config.path_graph_stats_fh = open(config.path_graph_stats_file, 'w')
    fields = ['locus', 'k', 'kmax', 'transfrags', 'kmers',
              'short_transfrags', 'lost_kmers', 'tot_expr', 'lost_expr',
              'lost_expr_frac', 'valid']
    print >>config.path_graph_stats_fh, '\t'.join(fields)

    # assembly gtf and bed files
    config.assembly_loss_gtf_fh = open(config.assembly_loss_gtf_file, 'w')
    config.assembly_gtf_fh = open(config.assembly_gtf_file, 'w')
    config.assembly_bed_fh = open(config.assembly_bed_file, 'w')

    # parse gtf file
    for interval, gtf_lines in GTF.parse_loci(open(config.transfrags_gtf_file)):
        chrom, start, end = interval
        logging.debug('Locus %s:%d-%d: ' % (chrom, start, end))
        assemble_locus(gtf_lines, config)

    # cleanup and close files
    config.assembly_gtf_fh.close()
    config.assembly_bed_fh.close()
    config.assembly_loss_gtf_fh.close()
    config.path_graph_stats_fh.close()
    config.splice_graph_gtf_fh.close()
    config.splice_bed_fh.close()
    Locus.close_bedgraphs(config.unresolved_bg_fhs)
    Locus.close_bedgraphs(config.resolved_bg_fhs)
