'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import bisect
import networkx as nx
import numpy as np

from base import Exon, Strand, TacoError, FLOAT_DTYPE
from locus import find_genomic_boundaries

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def aggregate_expression_data(transfrags, start, end):
    # add transfrag expression
    length = end - start
    expr_data = np.zeros(length, dtype=FLOAT_DTYPE)
    for t in transfrags:
        if t.is_ref:
            continue
        for exon in t.exons:
            astart = exon.start - start
            aend = exon.end - start
            expr_data[astart:aend] += t.expr
    return expr_data


def find_zero_sites(a, start=0):
    '''
    input: 'a' - numpy array of expression data
    '''
    # TODO: port to cython
    if a.shape[0] == 0:
        return []
    if a.shape[0] == 1:
        return []
    assert a[0] != 0
    assert a[-1] != 0
    boundaries = []
    for i in xrange(1, a.shape[0]):
        if (a[i-1] > 0) and (a[i] == 0):
            boundaries.append(start + i)
        elif (a[i-1] == 0) and (a[i] > 0):
            boundaries.append(start + i)
    return boundaries


def find_splice_sites(transfrags):
    boundaries = set()
    for t in transfrags:
        # add intron boundaries
        for start, end in t.iterintrons():
            boundaries.add(start)
            boundaries.add(end)
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


class SpliceGraph(object):

    def __init__(self):
        self.transfrags = None
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.boundaries = None
        self.expr_data = None
        self.sources = None
        self.sinks = None

    @staticmethod
    def create(transfrags):
        self = SpliceGraph()
        self.transfrags = transfrags

        chrom, start, end, strands = find_genomic_boundaries(transfrags)
        self.chrom = chrom
        self.start = start
        self.end = end
        assert len(strands) == 1
        self.strand = strands.pop()
        # aggregate transfrag expression
        expr_data = aggregate_expression_data(transfrags, start, end)
        self.expr_data = expr_data
        # define graph node boundaries
        boundaries = set([start, end])
        # nodes bounded by internal regions with zero expression
        boundaries.update(find_zero_sites(expr_data, start=start))
        # nodes bounded by introns
        boundaries.update(find_splice_sites(transfrags))
        self.boundaries = sorted(boundaries)
        return self

    def construct_graph(self):
        def add_node(G, n):
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=None)

        G = nx.DiGraph()
        for t in self.transfrags:
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [Exon(a, b) for (a, b) in
                     split_transfrag(t, self.boundaries)]
            if self.strand == Strand.NEG:
                nodes.reverse()
            # add nodes/edges to graph
            u = nodes[0]
            add_node(G, u)
            for i in xrange(1, len(nodes)):
                v = nodes[i]
                add_node(G, v)
                G.add_edge(u, v)
                u = v
        # set graph attributes
        G.graph['boundaries'] = self.boundaries
        G.graph['expr_data'] = self.expr_data
        G.graph['chrom'] = self.chrom
        G.graph['start'] = self.start
        G.graph['end'] = self.end
        G.graph['strand'] = self.strand
        return G

    def split_connected_components(self):
        # get connected components of graph which represent independent genes
        # unconnected components are considered different genes
        Gsubs = nx.weakly_connected_component_subgraphs(self.G)
        print Gsubs
