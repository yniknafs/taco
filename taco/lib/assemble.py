'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import bisect
import logging
import networkx as nx
import numpy as np

from gtf import GTF
from base import Strand, TacoError
from bedgraph import array_to_bedgraph
from transfrag import Transfrag

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
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
            self.strand_transfrags[t.strand].append(t)
            self._add_transfrag(t)
        return self

    def _add_transfrag(self, t):
        if t.is_ref:
            # TODO: handle ref transcripts
            return
        for exon in t.exons:
            start = exon.start - self.boundaries[0]
            end = exon.end - self.boundaries[0]
            self.expr_data[t.strand, start:end] += t.expr
            self.strand_data[t.strand, start:end] = True

    def _check_strand_ambiguous(self, t):
        '''
        Checks list of nodes for strandedness. If strand is unambiguous,
        then return pos or neg strand. If ambiguous, return unstranded.
        '''
        strands = [False, False]
        for exon in t.exons:
            start = exon.start - self.boundaries[0]
            end = exon.end - self.boundaries[0]
            for s in (Strand.POS, Strand.NEG):
                if self.strand_data[s, start:end].any():
                    strands[s] = True
                if self.expr_data[s, start:end].sum() > 0:
                    strands[s] = True
        if strands[Strand.POS] and strands[Strand.NEG]:
            return Strand.NA
        elif strands[Strand.POS]:
            return Strand.POS
        elif strands[Strand.NEG]:
            return Strand.NEG
        return Strand.NA

    def impute_unknown_strands(self):
        # iteratively predict strand of unstranded transfrags
        # stop when no new transfrag strands can be imputed
        iterations = 0
        num_resolved = 0
        while(len(self.strand_transfrags[Strand.NA]) > 0):
            resolved = []
            unresolved = []
            for t in self.strand_transfrags[Strand.NA]:
                new_strand = self._check_strand_ambiguous(t)
                if new_strand != Strand.NA:
                    resolved.append((t, new_strand))
                    num_resolved += 1
                else:
                    unresolved.append(t)
            # break when no new transfrags could have strand predicted
            if len(resolved) == 0:
                break
            # add resolved transfrags to new strand
            for t, new_strand in resolved:
                t.strand = new_strand
                self._add_transfrag(t)
                self.strand_transfrags[t.strand].append(t)
            self.strand_transfrags[Strand.NA] = unresolved
            iterations += 1

        # clear unstranded data and re-add unresolved transcripts
        if num_resolved > 0:
            self.strand_data[Strand.NA, :] = False
            self.expr_data[Strand.NA, :] = 0.0
            for t in self.strand_transfrags[Strand.NA]:
                self._add_transfrag(t)

        if iterations > 0:
            logging.debug('predict_unknown_strands: %d iterations' %
                          iterations)
        return num_resolved

    def get_expr_data(self, start, end, strand):
        if ((start < self.boundaries[0]) or (end > self.boundaries[-1])):
            raise TacoError('start/end are out of bounds')
        start = start - self.boundaries[0]
        end = end - self.boundaries[0]
        return self.expr_data[strand, start:end]

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
        boundaries = find_boundaries(transfrags)
        G = nx.DiGraph()

        # add transcripts
        for t in transfrags:
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [n for n in split_transfrag(t, boundaries)]
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
    def open_bedgraphs(file_prefix):
        bgfilehs = []
        for s in [Strand.POS, Strand.NEG, Strand.NA]:
            filename = '%s.%s.bedgraph' % (file_prefix, Strand.NAMES[s])
            bgfilehs.append(open(filename, 'w'))
        return bgfilehs

    @staticmethod
    def close_bedgraphs(bgfilehs):
        for fileh in bgfilehs:
            fileh.close()

    def get_bedgraph_data(self):
        '''
        Returns data in tuples
            chrom, start, end, strand, expression
        '''
        for strand in (Strand.POS, Strand.NEG, Strand.NA):
            for start, end, expr in array_to_bedgraph(self.expr_data[strand]):
                yield (self.chrom, start + self.boundaries[0],
                       end + self.boundaries[0], strand, expr)

    def write_bedgraph(self, bgfilehs):
        '''
        bgfiledict: dictionary structure containing file handles opened
                    for writing obtained using Locus.open_bedgraph()
        '''
        for tup in self.get_bedgraph_data():
            chrom, start, end, strand, expr = tup
            if expr > 0:
                line = '\t'.join(map(str, [chrom, start, end, expr]))
                print >>bgfilehs[strand], line


def assemble(gtf_file, output_dir):
    # setup bedgraph output files
    file_prefix = os.path.join(output_dir, 'loci.unresolved')
    raw_bgfilehd = Locus.open_bedgraphs(file_prefix)
    file_prefix = os.path.join(output_dir, 'loci.resolved')
    resolved_bgfilehd = Locus.open_bedgraphs(file_prefix)

    # parse gtf file
    for interval, gtf_lines in GTF.parse_loci(open(gtf_file)):
        chrom, start, end = interval
        print chrom, start, end, len(gtf_lines)
        t_dict = Transfrag.parse_gtf(gtf_lines)
        locus = Locus.create(t_dict.values())
        logging.debug('Locus %s:%d-%d: '
                      '%d transfrags (+: %d, -: %d, .: %d)' %
                      (chrom, start, end, len(t_dict),
                       len(locus.strand_transfrags[Strand.POS]),
                       len(locus.strand_transfrags[Strand.NEG]),
                       len(locus.strand_transfrags[Strand.NA])))
        # write raw bedgraph files
        locus.write_bedgraph(raw_bgfilehd)
        # resolve unstranded transcripts
        num_resolved = locus.impute_unknown_strands()
        if num_resolved > 0:
            logging.debug('Locus %s:%d-%d: %d '
                          'resolved (+: %d, -: %d, .: %d)' %
                          (chrom, start, end, num_resolved,
                           len(locus.strand_transfrags[Strand.POS]),
                           len(locus.strand_transfrags[Strand.NEG]),
                           len(locus.strand_transfrags[Strand.NA])))
        # write bedgraph files after strand resolved
        locus.write_bedgraph(resolved_bgfilehd)

    # close bedgraph files
    Locus.close_bedgraphs(raw_bgfilehd)
    Locus.close_bedgraphs(resolved_bgfilehd)
