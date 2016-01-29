'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import collections
from multiprocessing import Process, JoinableQueue, Value, Lock

from gtf import GTF
from base import Strand
from transfrag import Transfrag
from locus import Locus
from path_graph import choose_k, create_path_graph, smooth_graph, \
    reconstruct_path
from path_finder import find_suboptimal_paths
from bx.cluster import ClusterTree


__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class LockValue(object):
    def __init__(self, initval=0):
        self.val = Value('L', initval)
        self.lock = Lock()

    def next(self):
        with self.lock:
            cur_val = self.val.value
            self.val.value += 1
            return cur_val


class Isoform(object):
    __slots__ = ('expr', 'path', 'gene_id', 'tss_id')

    def __init__(self, path, expr):
        self.expr = expr
        self.path = path
        self.gene_id = -1
        self.tss_id = -1


class Config(object):
    def __init__(self, **kwargs):
        self.gtf_file = None
        self.unresolved_bg_files = None
        self.resolved_bg_files = None
        self.expr_h5_file = None
        self.chrom_sizes_file = None
        self.node_gtf_file = None
        self.assembly_gtf_file = None
        self.assembly_bed_file = None
        self.min_path_length = 400
        self.frac_isoform = 0.01
        self.max_isoforms = 100
        self.guided_strand = False
        self.guided_ends = False
        self.guided_assembly = False
        for k, v in kwargs.iteritems():
            setattr(self, k, v)


def get_gtf_features(chrom, strand, exons, locus_id, gene_id, tss_id,
                     transcript_id, expr, frac):
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
    f.score = int(round(1000.0 * frac))
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'expr': '%.3f' % expr,
               'frac': '%.3f' % frac}
    f.attrs.update(attr_dict)
    yield f
    for e in exons:
        f = GTF.Feature()
        f.seqid = chrom
        f.source = 'taco'
        f.feature = 'exon'
        f.start = e.start
        f.end = e.end
        f.score = int(round(1000.0 * frac))
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


def annotate_gene_and_tss_ids(isoforms, strand,
                              gene_id_value_obj,
                              tss_id_value_obj):
    # cluster paths to determine gene ids
    cluster_tree = ClusterTree(0, 1)
    # map tss positions to unique ids
    tss_pos_id_map = {}
    for i, isoform in enumerate(isoforms):
        start = isoform.path[0].start
        end = isoform.path[-1].end
        # cluster transcript coordinates
        cluster_tree.insert(start, end, i)
        # map TSS positions to IDs
        tss_pos = end if strand == Strand.NEG else start
        if tss_pos not in tss_pos_id_map:
            tss_id = tss_id_value_obj.next()
            tss_pos_id_map[tss_pos] = tss_id
        else:
            tss_id = tss_pos_id_map[tss_pos]
        isoform.tss_id = tss_id
    # retrieve transcript clusters and assign gene ids
    for start, end, indexes in cluster_tree.getregions():
        gene_id = gene_id_value_obj.next()
        for i in indexes:
            isoforms[i].gene_id = gene_id


def assemble_isoforms(sgraph, min_path_length, frac_isoform, max_isoforms):
    # choose a k-mer size 'k' based on a desired minimum path length
    k = choose_k(sgraph.transfrags,
                 sgraph.node_bounds,
                 min_path_length=min_path_length)
    K, lost_paths = create_path_graph(sgraph, k)
    # smooth kmer graph
    smooth_graph(K)
    # find up to 'max_isoforms' paths through graph
    logging.debug('%s:%d-%d[%s] finding paths in k=%d graph (%d nodes)' %
                  (sgraph.chrom, sgraph.start, sgraph.end,
                   Strand.to_gtf(sgraph.strand), k, len(K)))
    isoforms = []
    id_kmer_map = K.graph['id_kmer_map']
    for kmer_path, expr in \
        find_suboptimal_paths(K, K.graph['source'], K.graph['sink'],
                              frac_isoform, max_isoforms):
        # reconstruct path
        path = reconstruct_path(kmer_path, id_kmer_map, sgraph.strand)
        logging.debug("\texpr=%f length=%d" % (expr, len(path)))
        # add to path list
        isoforms.append(Isoform(path, expr))
    return isoforms


def assemble_gene(sgraph, locus_id_str, config):
    # run isoform path finding algorithm
    isoforms = assemble_isoforms(sgraph, config.min_path_length,
                                 config.frac_isoform, config.max_isoforms)
    # determine gene ids and tss ids
    annotate_gene_and_tss_ids(isoforms, sgraph.strand,
                              config.gene_id_value_obj,
                              config.tss_id_value_obj)

    # bin transcripts by gene id
    gene_isoform_dict = collections.defaultdict(lambda: [])
    for isoform in isoforms:
        gene_isoform_dict[isoform.gene_id].append(isoform)

    # write output
    for gene_isoforms in gene_isoform_dict.itervalues():
        # compute total expression
        total_expr = sum(isoform.expr for isoform in gene_isoforms)
        total_expr = max(1e-8, total_expr)
        # create GTF features for each transcript path
        for isoform in gene_isoforms:
            # assign transcript id
            t_id = config.t_id_value_obj.next()
            # get strings for each id
            t_id_str = "TU%d" % t_id
            tss_id_str = "TSS%d" % (isoform.tss_id)
            gene_id_str = "G%d" % (isoform.gene_id)
            # compute isoform fractions
            frac = isoform.expr / total_expr
            # write to GTF
            for f in get_gtf_features(chrom=sgraph.chrom,
                                      strand=sgraph.strand,
                                      exons=isoform.path,
                                      locus_id=locus_id_str,
                                      gene_id=gene_id_str,
                                      tss_id=tss_id_str,
                                      transcript_id=t_id_str,
                                      expr=isoform.expr,
                                      frac=frac):
                print >>config.assembly_gtf_fh, str(f)
            # write to BED
            name = "%s|%s(%.1f)" % (gene_id_str, t_id_str, isoform.expr)
            fields = write_bed(sgraph.chrom, name, sgraph.strand,
                               int(round(1000.0*frac)), isoform.path)
            print >>config.bed_fh, '\t'.join(fields)


def assemble_locus(gtf_lines, config):
    t_dict = Transfrag.parse_gtf(gtf_lines)
    locus = Locus.create(t_dict.values(),
                         config.guided_strand,
                         config.guided_ends,
                         config.guided_assembly)
    logging.debug('Locus %s:%d-%d: '
                  '%d transfrags (+: %d, -: %d, .: %d)' %
                  (locus.chrom, locus.start, locus.end, len(t_dict),
                   len(locus.get_transfrags(Strand.POS)),
                   len(locus.get_transfrags(Strand.NEG)),
                   len(locus.get_transfrags(Strand.NA))))
    # write raw bedgraph files
    locus.write_bedgraph(config.unresolved_bg_fhs)
    # resolve unstranded transcripts
    num_resolved = locus.impute_unknown_strands()
    if num_resolved > 0:
        logging.debug('Locus %s:%d-%d: %d '
                      'resolved (+: %d, -: %d, .: %d)' %
                      (locus.chrom, locus.start, locus.end, num_resolved,
                       len(locus.get_transfrags(Strand.POS)),
                       len(locus.get_transfrags(Strand.NEG)),
                       len(locus.get_transfrags(Strand.NA))))
    # write bedgraph files after strand resolved
    locus.write_bedgraph(config.resolved_bg_fhs)
    # write splice junctions
    locus.write_splice_bed(config.splice_bed_fh)
    # write expression array
    locus.write_expression_hdf5(config.expr_h5fh)
    # convert to stranded locus objects
    locus_id_str = "L%d" % (config.locus_id_value_obj.next())
    for sgraph in locus.create_splice_graphs():
        for f in sgraph.get_node_gtf():
            print >>config.node_gtf_fh, str(f)
        assemble_gene(sgraph, locus_id_str, config)


def assemble(**kwargs):
    config = Config(**kwargs)
    # setup bedgraph output files
    config.unresolved_bg_fhs = []
    for s, filename in config.unresolved_bg_files:
        config.unresolved_bg_fhs.append(open(filename, 'w'))
    config.resolved_bg_fhs = []
    for s, filename in config.resolved_bg_files:
        config.resolved_bg_fhs.append(open(filename, 'w'))
    # setup junction bed file
    config.splice_bed_fh = Locus.open_splice_bed(config.splice_bed_file)
    # setup expression hdf5
    config.expr_h5fh = Locus.open_expression_hdf5(config.expr_h5_file,
                                                  config.chrom_sizes_file)
    # setup node gtf file
    config.node_gtf_fh = open(config.node_gtf_file, 'w')
    # assembly gtf and bed files
    config.assembly_gtf_fh = open(config.assembly_gtf_file, 'w')
    config.bed_fh = open(config.assembly_bed_file, 'w')
    # shared memory values (for parallelism)
    config.locus_id_value_obj = LockValue(1)
    config.gene_id_value_obj = LockValue(1)
    config.tss_id_value_obj = LockValue(1)
    config.t_id_value_obj = LockValue(1)

    # parse gtf file
    for interval, gtf_lines in GTF.parse_loci(open(config.gtf_file)):
        chrom, start, end = interval
        logging.debug('Locus %s:%d-%d: ' % (chrom, start, end))
        assemble_locus(gtf_lines, config)

    # close assembly files
    config.assembly_gtf_fh.close()
    config.bed_fh.close()
    # close node file
    config.node_gtf_fh.close()
    # close expression
    config.expr_h5fh.close()
    # close splice junction bed file
    config.splice_bed_fh.close()
    # close bedgraph files
    Locus.close_bedgraphs(config.unresolved_bg_fhs)
    Locus.close_bedgraphs(config.resolved_bg_fhs)
