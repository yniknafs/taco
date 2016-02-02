'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import collections
import bisect
from itertools import chain
import networkx as nx
import numpy as np

from base import Exon, Strand, TacoError
from gtf import GTF
from dtypes import FLOAT_DTYPE
from cChangePoint import find_threshold_points
from changepoint import run_changepoint
from bx.intersection import Interval, IntervalTree

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


# graph node attributes
class Node:
    IS_START = 'is_start'
    IS_STOP = 'is_stop'
    IS_REF = 'is_ref'

    @staticmethod
    def init():
        return {Node.IS_START: False,
                Node.IS_STOP: False,
                Node.IS_REF: False}

    @staticmethod
    def add(G, n, is_ref=False):
        if n not in G:
            G.add_node(n, attr_dict=Node.init())
        if is_ref:
            G.node[n][Node.IS_REF] = True


def _splice_graph_add_transfrag(G, t, node_bounds):
    # split exons that cross boundaries and get the
    # nodes that made up the transfrag
    nodes = [Exon(a, b) for (a, b) in
             split_transfrag(t, node_bounds)]
    if t.strand == Strand.NEG:
        nodes.reverse()
    # add nodes/edges to graph
    u = nodes[0]
    Node.add(G, u, is_ref=t.is_ref)
    for i in xrange(1, len(nodes)):
        v = nodes[i]
        Node.add(G, v, is_ref=t.is_ref)
        G.add_edge(u, v)
        u = v
    return G


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


def split_transfrag(t, boundaries):
    '''
    output: (generator) tuples (a,b) reflecting nodes
    '''
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
        self.guided_ends = False
        self.guided_assembly = False
        self.transfrags = []
        self.ref_transfrags = []
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.tree_starts = IntervalTree()
        self.tree_ends = IntervalTree()
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        self.start_sites = set()
        self.stop_sites = set()
        self.expr_data = None
        self.node_bounds = None
        self.G = None

    @staticmethod
    def create(transfrags,
               guided_ends=False,
               guided_assembly=False):
        self = SpliceGraph()
        self.guided_ends = guided_ends
        self.guided_assembly = guided_assembly
        self._add_transfrags(transfrags)
        self.expr_data = self._compute_expression()
        self.node_bounds = self._find_node_boundaries()
        self.G = self._create_splice_graph()
        self._mark_start_stop_nodes()
        return self

    def recreate(self):
        self._add_transfrags(chain(self.transfrags, self.ref_transfrags))
        self.expr_data = self._compute_expression()
        self.node_bounds = self._find_node_boundaries()
        self.G = self._create_splice_graph()
        self._mark_start_stop_nodes()

    def _add_transfrags(self, transfrags):
        self.transfrags = []
        self.ref_transfrags = []
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        self.tree_starts = IntervalTree()
        self.tree_ends = IntervalTree()
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
            if t.is_ref:
                self.ref_start_sites.add(t.txstart)
                self.ref_stop_sites.add(t.txstop)
                self.ref_transfrags.append(t)
            else:
                self._add_to_interval_trees(t)
                self.transfrags.append(t)
        self.ref_start_sites = sorted(self.ref_start_sites)
        self.ref_stop_sites = sorted(self.ref_stop_sites)

    def _compute_expression(self):
        expr_data = np.zeros(self.end - self.start, dtype=FLOAT_DTYPE)
        for t in self.transfrags:
            for exon in t.exons:
                astart = exon.start - self.start
                aend = exon.end - self.start
                expr_data[astart:aend] += t.expr
        return expr_data

    def _find_node_boundaries(self):
        node_bounds = set((self.start, self.end))
        node_bounds.update(self.start_sites)
        node_bounds.update(self.stop_sites)
        # nodes bounded by regions where expression changes to/from zero
        node_bounds.update(find_threshold_points(self.expr_data, self.start))
        # nodes bounded by introns
        for t in self.transfrags:
            node_bounds.update(t.itersplices())
        if self.guided_ends or self.guided_assembly:
            for t in self.ref_transfrags:
                if self.guided_ends:
                    node_bounds.update((t.start, t.end))
                if self.guided_assembly:
                    node_bounds.update(t.itersplices())
        return sorted(node_bounds)

    def _create_splice_graph(self):
        '''returns networkx DiGraph object'''
        G = nx.DiGraph()
        for t in self.itertransfrags():
            _splice_graph_add_transfrag(G, t, self.node_bounds)
        return G

    def _mark_start_stop_nodes(self):
        G = self.G
        # get all leaf nodes
        for n, nd in G.nodes_iter(data=True):
            if G.in_degree(n) == 0:
                nd[Node.IS_START] = True
            if G.out_degree(n) == 0:
                nd[Node.IS_STOP] = True
        # mark change points
        change_points = set()
        change_points.update(set((Node.IS_START, x) for x in self.start_sites))
        change_points.update(set((Node.IS_STOP, x) for x in self.stop_sites))
        if self.guided_ends:
            change_points.update((Node.IS_START, x)
                                 for x in self.ref_start_sites)
            change_points.update((Node.IS_STOP, x)
                                 for x in self.ref_stop_sites)
        node_bounds = self.node_bounds
        strand = self.strand
        for direction, pos in change_points:
            if ((direction == Node.IS_STOP and strand == Strand.NEG) or
                (direction == Node.IS_START and strand != Strand.NEG)):
                bisect_func = bisect.bisect_right
            else:
                bisect_func = bisect.bisect_left
            i = bisect_func(node_bounds, pos)
            n = (node_bounds[i-1], node_bounds[i])
            assert n in G
            G.node[n][direction] = True

    def _add_to_interval_trees(self, t):
        istart = Interval(t.start, t.start+1, value=t)
        iend = Interval(t.end-1, t.end, value=t)
        self.tree_starts.insert_interval(istart)
        self.tree_ends.insert_interval(iend)

    def _trim_change_point(self, cp_pos, cp_start, cp_end, cp_sign):
        # search for matches in change point interval
        num_trimmed = 0
        if cp_sign < 0:
            hits = self.tree_ends.find(cp_pos, cp_end)
            for hit in hits:
                t = hit.value
                # last exon start cannot overlap interval
                last_exon_start = t.exons[-1][0]
                if cp_pos <= last_exon_start <= cp_end:
                    continue
                # trim end left
                t.exons[-1] = Exon(last_exon_start, cp_pos)
                num_trimmed += 1
        else:
            hits = self.tree_starts.find(cp_start, cp_pos)
            for hit in hits:
                t = hit.value
                # first exon end cannot overlap interval
                first_exon_end = t.exons[0][1]
                if cp_start <= first_exon_end <= cp_pos:
                    continue
                # trim start right
                t.exons[0] = Exon(cp_pos, first_exon_end)
                num_trimmed += 1
        logging.debug('\tfound %d hits and trimmed %d' % (len(hits), num_trimmed))

    def detect_change_points(self, trim=True, *args, **kwargs):
        '''
        trim: if true, will trim transfrags around change points
        *args, **kwargs: passed directly to 'run_changepoint'
        '''
        for node in self.G:
            expr_data = self.get_expr_data(node.start, node.end)
            changepts = run_changepoint(expr_data, *args, **kwargs)
            # add start and stop sites
            for cp in changepts:
                cp_left = cp.index - cp.dist_left
                cp_right = cp.index + cp.dist_right
                cp_pos = node.start + cp.index
                cp_start = node.start + cp_left
                cp_end = node.start + cp_right
                logging.debug('\tchange point: %s:%d(%d-%d)[%s] node(%d-%d) '
                              'sign=%f pval=%f fc=%f' %
                              (self.chrom, cp_pos, cp_start, cp_end,
                               Strand.to_gtf(self.strand),
                               node.start, node.end,
                               cp.sign, cp.pvalue, cp.foldchange))
                if ((self.strand != Strand.NEG and cp.sign < 0) or
                    (self.strand == Strand.NEG and cp.sign > 0)):
                    self.stop_sites.add(cp_pos)
                else:
                    self.start_sites.add(cp_pos)
                if trim:
                    self._trim_change_point(cp_pos, cp_start, cp_end, cp.sign)

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
        for t in self.itertransfrags():
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

    def itertransfrags(self):
        if self.guided_assembly:
            return chain(self.transfrags, self.ref_transfrags)
        else:
            return iter(self.transfrags)

    def get_start_stop_nodes(self):
        start_nodes = set()
        stop_nodes = set()
        for n, nd in self.G.nodes_iter(data=True):
            if nd[Node.IS_START]:
                start_nodes.add(n)
            if nd[Node.IS_STOP]:
                stop_nodes.add(n)
        return start_nodes, stop_nodes

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
        for n in self.G:
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
