'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import bisect
import networkx as nx
import numpy as np

from base import Exon, Strand, TacoError
from dtypes import FLOAT_DTYPE
from cChangePoint import find_change_points

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


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



def aggregate_expression_data(t, expr_data, start):
    # add transfrag expression
    for exon in t.exons:
        astart = exon.start - start
        aend = exon.end - start
        expr_data[astart:aend] += t.expr
    return expr_data


def find_splice_sites(t):
    a = []
    for start, end in t.iterintrons():
        a.append(start)
        a.append(end)
    return a


def get_start_stop_sites(t):
    t_start = t.exons[0][0]
    t_stop = t.exons[-1][1]
    if t.strand == Strand.NEG:
        t_start, t_stop = t_stop, t_start
    return t_start, t_stop




class SpliceGraph(object):

    def __init__(self):
        self.guided_ends = False
        self.guided_assembly = False
        self.transfrags = None
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.expr_data = None
        self.boundaries = None
        self.ref_start_sites = None
        self.ref_stop_sites = None
        self.start_sites = None
        self.stop_sites = None

    @staticmethod
    def create(transfrags, guided_ends=False, guided_assembly=False):
        self = SpliceGraph()
        self.guided_ends = guided_ends
        self.guided_assembly = guided_assembly
        self.transfrags = transfrags
        # first pass through transfrags to define graph boundaries
        for t in transfrags:
            if self.chrom is None:
                self.chrom = t.chrom
            elif self.chrom != t.chrom:
                raise TacoError('chrom mismatch')
            if self.strand is None:
                self.strand = t.strand
            elif self.strand != t.strand:
                raise TacoError('strand mismatch')
            if self.start is None:
                self.start = t.start
            else:
                self.start = min(t.start, self.start)
            if self.end is None:
                self.end = t.end
            else:
                self.end = max(t.end, self.end)
        # second pass through transfrags to aggregate expression and
        # determine node boundaries
        length = self.end - self.start
        self.expr_data = np.zeros(length, dtype=FLOAT_DTYPE)
        self.boundaries = set([self.start, self.end])
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        for t in transfrags:
            if t.is_ref:
                t_start, t_stop = get_start_stop_sites(t)
                self.ref_start_sites.add(t_start)
                self.ref_stop_sites.add(t_stop)
                if guided_assembly:
                    self.boundaries.update(find_splice_sites(t))
                if guided_ends:
                    self.boundaries.update((t_start, t_stop))
            else:
                aggregate_expression_data(t, self.expr_data, self.start)
                self.boundaries.update(find_splice_sites(t))
        # nodes bounded by internal regions with zero expression
        zero_sites = find_change_points(self.expr_data, start=self.start)
        self.boundaries.update(zero_sites)
        self.boundaries = sorted(self.boundaries)
        return self

    def get_change_point_data(self):
        def get_region(start, end):
            # subset expression array
            astart = start - self.start
            aend = end - self.start
            expr_data = self.expr_data[astart:aend]
            # extract reference start/stop sites that occur within the region
            ref_starts = []
            for x in self.ref_start_sites:
                if (x >= start) and (x < end):
                    ref_starts.append(x - start)
            ref_starts.sort()
            ref_stops = []
            for x in self.ref_stop_sites:
                if (x > start) and (x <= end):
                    ref_stops.append(x - start)
            ref_stops.sort()
            return expr_data, ref_starts, ref_stops

        assert len(self.boundaries) >= 2
        bstart = self.boundaries[0]
        for bend in self.boundaries[1:]:
            expr_data, ref_starts, ref_stops = \
                get_region(bstart, bend)
            if expr_data.sum() == 0:
                continue
            line = '\t'.join([self.chrom, str(bstart), str(bend),
                              Strand.to_gtf(self.strand),
                              ','.join(map(str, ref_starts)),
                              ','.join(map(str, ref_stops)),
                              ','.join(map(str, expr_data))])
            yield line
            bstart = bend

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
