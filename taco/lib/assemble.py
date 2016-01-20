'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import collections
import bisect
import logging
import networkx as nx
import numpy as np

from gtf import GTF
from base import Strand, TacoError
from transfrag import Transfrag

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def find_boundaries(transfrags):
    '''
    input: list of transfrags
    output: sorted list of node boundaries
    '''
    boundaries = set()
    minstart = None
    maxend = None
    # first add introns to the graph and keep track of
    # all intron boundaries
    for t in transfrags:
        minstart = t.start if minstart is None else min(t.start, minstart)
        maxend = t.end if maxend is None else max(t.end, maxend)
        # add intron boundaries
        for start, end in t.iterintrons():
            boundaries.add(start)
            boundaries.add(end)
    # add the locus boundaries
    boundaries.add(minstart)
    boundaries.add(maxend)
    # sort the intron boundary positions
    return sorted(boundaries)


def split_transfrag(t, boundaries):
    '''
    input: transfrag
    output: (generator) tuples (a,b) reflecting transfrag nodes
    '''
    if t.start == t.end:
        return
    for exon in t.exons:
        # find the indexes into the splice sites list that border the exon.
        start_ind = bisect.bisect_right(boundaries, exon.start)
        end_ind = bisect.bisect_left(boundaries, exon.end)
        if start_ind == end_ind:
            yield boundaries[start_ind-1], boundaries[start_ind]
        else:
            # all the splice sites in between the exon borders must overlap
            for j in xrange(start_ind-1, end_ind):
                yield boundaries[j], boundaries[j+1]


class Locus(object):

    def __init__(self, guided=False):
        self.chrom = None
        self.boundaries = None
        self.expr_data = None
        self.strand_data = None
        self.sources = set()
        self.sinks = set()
        self.strand_transfrags = [[], [], []]
        self.guided = guided

    @staticmethod
    def create(transfrags, guided=False):
        self = Locus(guided=guided)
        self.boundaries = find_boundaries(transfrags)
        length = self.boundaries[-1] - self.boundaries[0]
        self.expr_data = np.zeros((3, length))
        self.strand_data = np.zeros((3, length), dtype=np.bool)

        for t in transfrags:
            if self.chrom is None:
                self.chrom = t.chrom
            elif self.chrom != t.chrom:
                logging.error('Locus.create: "chrom" mismatch')
                raise TacoError('Locus.create: transfrag chromosomes do '
                                'not match')
            self._add_transfrag(t)
            self.strand_transfrags[t.strand].append(t)
        return self

    def _add_transfrag(self, t):
        for n in split_transfrag(t, self.boundaries):
            start = n[0] - self.boundaries[0]
            end = n[1] - self.boundaries[0]
            if t.is_ref:
                # TODO: add starts/stops
                pass
            else:
                self.expr_data[t.strand, start:end] += t.expr
            self.strand_data[t.strand, start:end] = True
            print 't %s node %s' % (t._id, str(n))
        print 'expr strand %s %s' % (Strand.to_gtf(t.strand), ','.join(map(str, self.expr_data[t.strand])))

    def _remove_transfrag(self, t):
        for n in split_transfrag(t, self.boundaries):
            start = n[0] - self.boundaries[0]
            end = n[1] - self.boundaries[0]

            self.strand_data[t.strand]

            nd = self.node_data[n]
            nd.strands[t.strand] = False
            nd.samples[t.strand].remove(t.sample_id)
            nd.exprs[t.strand] -= t.expr

    def _check_strand_ambiguous(self, nodes):
        '''
        Checks list of nodes for strandedness. If strand is unambiguous,
        then return pos or neg strand. If ambiguous, return unstranded.
        '''
        strands = {GTF.POS_STRAND: False,
                   GTF.NEG_STRAND: False}
        for n in nodes:
            nd = self.node_data[n]
            if nd.strands['+'] or nd.exprs['+'] > 0:
                strands['+'] = True
            if nd.strands['-'] or nd.exprs['-'] > 0:
                strands['-'] = True
        if strands['+'] and strands['-']:
            return '.'
        elif strands['+']:
            return '+'
        elif strands['-']:
            return '-'
        return '.'

    def impute_unknown_strands(self):
        # iteratively predict strand of unstranded transfrags
        # stop when no new transfrag strands can be imputed
        iterations = 0
        num_resolved = 0
        while(len(self.strand_transfrags['.']) > 0):
            resolved = []
            unresolved = []
            for t in self.strand_transfrags['.']:
                nodes = list(split_exons(t, self.boundaries))
                new_strand = self._check_strand_ambiguous(nodes)
                if new_strand != '.':
                    resolved.append((t, new_strand))
                    num_resolved += 1
                else:
                    unresolved.append(t)
            # break when no new transfrags could have strand predicted
            if len(resolved) == 0:
                break
            # clear Locus data for transfrags with unknown strand and re-add
            # them to the Locus
            for t, new_strand in resolved:
                self._remove_transfrag(t)
                t.strand = new_strand
                self._add_transfrag(t)
                self.strand_transfrags[t.strand].append(t)
            self.strand_transfrags['.'] = unresolved
            iterations += 1

        if iterations > 0:
            logging.debug('predict_unknown_strands: %d iterations' %
                          iterations)
        return num_resolved

    def create_directed_graph(self, strand):
        '''
        build strand-specific graph
        '''
        def add_node(G, n, expr):
            """add node to graph"""
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=0.0)
            nd = G.node[n]
            nd['expr'] += expr

        # initialize transcript graph
        transfrags = self.strand_transfrags[strand]
        boundaries = find_exon_boundaries(transfrags)
        G = nx.DiGraph()

        # add transcripts
        for t in transfrags:
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [n for n in split_exons(t, boundaries)]
            if strand == '-':
                nodes.reverse()
            # add nodes/edges to graph
            u = nodes[0]
            add_node(G, u, t.expr)
            for i in xrange(1, len(nodes)):
                v = nodes[i]
                add_node(G, v, t.expr)
                G.add_edge(u, v)
                u = v

        # set graph attributes
        G.graph['boundaries'] = boundaries
        G.graph['strand'] = strand
        return G

    def create_directed_graphs(self):
        for strand, transfrags in self.strand_transfrags.iteritems():
            # create strand-specific directed graph
            G = self.create_directed_graph(strand)
            # collapse consecutive nodes in graph
            H, node_chain_map = collapse_strand_specific_graph(G, introns=True)
            # get connected components of graph which represent independent genes
            # unconnected components are considered different genes
            Gsubs = nx.weakly_connected_component_subgraphs(H)

    @staticmethod
    def open_bedgraph(file_prefix):
        attrs = ('expression', 'recurrence')
        strand_names = {'+': 'pos', '-': 'neg', '.': 'none'}
        filehs = {}
        for a in attrs:
            filehs[a] = {}
            for s in ('+', '-', '.'):
                filename = '%s.%s.%s.bedgraph' % (file_prefix, a,
                                                  strand_names[s])
                filehs[a][s] = open(filename, 'w')
        return filehs

    @staticmethod
    def close_bedgraph(filehs):
        for adict in filehs.itervalues():
            for fileh in adict.itervalues():
                fileh.close()

    def get_bedgraph_data(self):
        '''
        Returns node attribute data in tuples
            chrom, start, end, strand, total expression, sample recurrence
        '''
        for n in sorted(self.node_data):
            nd = self.node_data[n]
            for strand in ('+', '-', '.'):
                yield (self.chrom, n[0], n[1], strand,
                       nd.exprs[strand], len(nd.samples[strand]))

    def write_bedgraph(self, bgfiledict):
        '''
        bgfiledict: dictionary structure containing file handles opened
                    for writing obtained using Locus.open_bedgraph()
        '''
        for tup in self.get_bedgraph_data():
            chrom, start, end, strand, expr, recur = tup
            if expr > 0:
                line = '\t'.join(map(str, [chrom, start, end, expr]))
                print >>bgfiledict['expression'][strand], line
            if recur > 0:
                line = '\t'.join(map(str, [chrom, start, end, recur]))
                print >>bgfiledict['recurrence'][strand], line


def impute_strand(a, r):
    # setup bedgraph output files
    # file_prefix = os.path.join(a.output_dir, 'loci.unresolved')
    # raw_bgfilehd = Locus.open_bedgraph(file_prefix)
    # file_prefix = os.path.join(a.output_dir, 'loci.resolved')
    # resolved_bgfilehd = Locus.open_bedgraph(file_prefix)
    import sys

    # parse gtf file
    ignore_ref = not a.guided
    for interval, gtf_lines in GTF.parse_loci(open(r.transfrags_gtf_file)):
        chrom, start, end = interval
        t_dict = Transfrag.parse_gtf(gtf_lines, ignore_ref)
        locus = Locus.create(t_dict.values())
        logging.debug('Locus %s:%d-%d: '
                      '%d transfrags (+: %d, -: %d, .: %d)' %
                      (chrom, start, end, len(t_dict),
                       len(locus.strand_transfrags['+']),
                       len(locus.strand_transfrags['-']),
                       len(locus.strand_transfrags['.'])))
        sys.exit(0)
        # write bedgraph files for expression/recurrence data
        locus.write_bedgraph(raw_bgfilehd)

        # resolve unstranded transcripts
        num_resolved = locus.impute_unknown_strands()
        if num_resolved > 0:
            logging.debug('Locus %s:%d-%d: %d '
                          'resolved (+: %d, -: %d, .: %d)' %
                          (chrom, start, end, num_resolved,
                           len(locus.strand_transfrags['+']),
                           len(locus.strand_transfrags['-']),
                           len(locus.strand_transfrags['.'])))

        # write bedgraph files after strand resolved
        locus.write_bedgraph(resolved_bgfilehd)

    # close bedgraph files
    Locus.close_bedgraph(raw_bgfilehd)
    Locus.close_bedgraph(resolved_bgfilehd)

    # update status and write to file
    # self.status.assemble = True
    # self.status.write(self.results.status_file)
    # return EXIT_SUCCESS
