'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import collections
import bisect
import networkx as nx
import numpy as np

from base import Exon, Strand, TacoError
from gtf import GTF
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


def _array_subset(a, start, end):
    '''
    return elements of sorted array 'a' where all elements of
    'start' < a[i] < 'end'
    '''
    assert start <= end
    if len(a) == 0:
        return []
    i = bisect.bisect_right(a, start)
    j = bisect.bisect_right(a, end)
    return a[i:j]


def find_node_boundaries(transfrags,
                         guided_ends=False,
                         guided_assembly=False):
    node_bounds = set()
    min_start = None
    max_end = None
    for t in transfrags:
        if t.is_ref:
            if guided_ends:
                node_bounds.update((t.start, t.end))
            if guided_assembly:
                node_bounds.update(t.itersplices())
        else:
            node_bounds.update(t.itersplices())
        min_start = t.start if min_start is None else min(t.start, min_start)
        max_end = t.end if max_end is None else max(t.end, max_end)
    node_bounds.add(min_start)
    node_bounds.add(max_end)
    return sorted(node_bounds)


def split_transfrag(t, boundaries):
    '''
    output: (generator) tuples (a,b) reflecting nodes
    '''
    for exon in t.exons:
        # find the indexes into the splice sites list that border the exon.
        start_ind = bisect.bisect_right(boundaries, exon.start)
        end_ind = bisect.bisect_left(boundaries, exon.end)

        start_ind, end_ind
        if start_ind == end_ind:
            yield boundaries[start_ind-1], boundaries[start_ind]
        else:
            # all the splice sites in between the exon borders must overlap
            for j in xrange(start_ind-1, end_ind):
                yield boundaries[j], boundaries[j+1]


class SpliceGraph(object):
    def __init__(self):
        self.guided_ends = False
        self.guided_assembly = False
        self.transfrags = []
        self.G = None
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.expr_data = None
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        self.start_sites = set()
        self.stop_sites = set()
        self.node_bounds = None
        self.G = None

    def _find_node_boundaries(self):
        node_bounds = set((self.start, self.end))
        # nodes bounded by regions where expression changes to/from zero
        zero_points = find_change_points(self.expr_data, start=self.start)
        node_bounds.update(zero_points)
        # nodes bounded by introns
        for t in self.transfrags:
            if t.is_ref:
                if self.guided_ends:
                    node_bounds.update((t.start, t.end))
                if self.guided_assembly:
                    node_bounds.update(t.itersplices())
            else:
                node_bounds.update(t.itersplices())
        return sorted(node_bounds)

    def _create_splice_graph(self):
        '''returns networkx DiGraph object'''
        def add_node(G, n):
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=None)

        G = nx.DiGraph()
        for t in self.transfrags:
            if t.is_ref and not self.guided_assembly:
                continue
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [Exon(a, b) for (a, b) in
                     split_transfrag(t, self.node_bounds)]
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
        return G

    @staticmethod
    def create(transfrags,
               guided_ends=False,
               guided_assembly=False):
        self = SpliceGraph()
        self.guided_ends = guided_ends
        self.guided_assembly = guided_assembly

        for t in transfrags:
            if self.chrom is None:
                self.chrom = t.chrom
            elif self.chrom != t.chrom:
                raise TacoError('chrom mismatch')
            if self.strand is None:
                self.strand = t.strand
            elif self.strand != t.strand:
                raise TacoError('chrom mismatch')
            if self.start is None:
                self.start = t.start
            else:
                self.start = min(t.start, self.start)
            if self.end is None:
                self.end = t.end
            else:
                self.end = max(t.end, self.end)
            if t.is_ref:
                self.ref_start_sites.add(t.txstart)
                self.ref_stop_sites.add(t.txstop)
            self.transfrags.append(t)

        self.ref_start_sites = sorted(self.ref_start_sites)
        self.ref_stop_sites = sorted(self.ref_stop_sites)
        self.expr_data = np.zeros(self.end - self.start, dtype=FLOAT_DTYPE)
        for t in self.transfrags:
            if not t.is_ref:
                for exon in t.exons:
                    astart = exon.start - self.start
                    aend = exon.end - self.start
                    self.expr_data[astart:aend] += t.expr

        self.node_bounds = self._find_node_boundaries()
        self.G = self._create_splice_graph()
        return self

    def split(self):
        '''splits into weakly connected component subgraphs'''
        # get connected components of graph which represent independent genes
        # unconnected components are considered different genes
        Gsubs = list(nx.weakly_connected_component_subgraphs(self.G))
        if len(Gsubs) == 1:
            yield self
            return
        # map nodes to components
        node_subgraph_map = {}
        subgraph_transfrag_map = collections.defaultdict(list)
        for i, Gsub in enumerate(Gsubs):
            for n in Gsub:
                node_subgraph_map[n] = i
        # assign transfrags to components
        for t in self.transfrags:
            if t.is_ref and not self.guided_assembly:
                continue
            for n in split_transfrag(t, self.node_bounds):
                n = Exon(*n)
                subgraph_id = node_subgraph_map[n]
                subgraph_transfrag_map[subgraph_id].append(t)
                break
        # create new graphs using the separate components
        for subgraph_transfrags in subgraph_transfrag_map.itervalues():
            yield SpliceGraph.create(subgraph_transfrags,
                                     guided_ends=self.guided_ends,
                                     guided_assembly=self.guided_assembly)

    def get_transfrags(self):
        for t in self.transfrags:
            if t.is_ref and not self.guided_assembly:
                continue
            yield t

    def get_expr_data(self, start, end):
        assert start >= self.start
        assert end <= self.end
        if ((start < self.start) or (end > self.end)):
            m = ('query %d-%d outside locus bounds %d-%d' %
                 (start, end, self.start, self.end))
            raise TacoError(m)
        astart = start - self.start
        aend = end - self.start
        return self.expr_data[astart:aend]

    def get_node_gtf(self):
        locus_id = ('L_%s_%d_%d_%s' %
                    (self.chrom, self.start, self.end,
                     Strand.to_gtf(self.strand)))
        # iterate through locus and return change point data
        for n in sorted(self.G):
            expr_data = self.get_expr_data(*n)
            ref_starts = _array_subset(self.ref_start_sites, *n)
            ref_stops = _array_subset(self.ref_stop_sites, *n)
            # return gtf feature for each node
            f = GTF.Feature()
            f.seqid = self.chrom
            f.source = 'taco'
            f.feature = 'node'
            f.start = n[0]
            f.end = n[1]
            f.score = 0
            f.strand = Strand.to_gtf(self.strand)
            f.phase = '.'
            f.attrs = {'locus_id': locus_id,
                       'expr_min': str(expr_data.min()),
                       'expr_max': str(expr_data.max()),
                       'expr_mean': str(expr_data.mean()),
                       'ref_starts': ','.join(map(str, ref_starts)),
                       'ref_stops': ','.join(map(str, ref_stops))}
            yield f
