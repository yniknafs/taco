'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import collections
import bisect
from itertools import chain
from array import array
import numpy as np

from base import Exon, Strand, TacoError
from graph import Graph
from gtf import GTF
from dtypes import FLOAT_DTYPE
from cchangepoint import find_threshold_points
from changepoint import run_changepoint
from bx.intersection import Interval, IntervalTree
import cbisect

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class SGraph(Graph):
    def __init__(self):
        super(SGraph, self).__init__()
        self.is_start = [True, True]
        self.is_stop = [True, True]
        self.is_ref = [False, False]

    def add_node(self, node):
        if node not in self.node_id_map:
            node_id = super(SGraph, self).add_node(node)
            self.is_start.append(False)
            self.is_stop.append(False)
            self.is_ref.append(False)
        else:
            node_id = self.node_id_map[node]
        return node_id

    def add_path(self, path, is_ref):
        node_ids = super(SGraph, self).add_path(path)
        for node_id in node_ids:
            if is_ref:
                self.is_ref[node_id] = True
        return node_ids

    def set_start_stop_nodes(self):
        # get all leaf nodes
        for i in self.node_ids_iter():
            if len(self.preds[i]) == 0:
                self.is_start[i] = True
            if len(self.succs[i]) == 0:
                self.is_stop[i] = True


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
    exons must be sorted in increasing order

    output: (generator) tuples (a,b) reflecting nodes
    '''
    end_ind = 0
    for exon in t.exons:
        # find the indexes into the splice sites list that border the exon.
        start_ind = cbisect.bisect_right(boundaries, exon.start, lo=end_ind)
        end_ind = cbisect.bisect_left(boundaries, exon.end, lo=start_ind)
        # start_ind = bisect.bisect_right(boundaries, exon.start)
        # end_ind = bisect.bisect_left(boundaries, exon.end)
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
        return array('i', sorted(node_bounds))

    def get_path(self, t):
        node_id_map = self.G.node_id_map
        nodes = [node_id_map[n] for n in split_transfrag(t, self.node_bounds)]
        if self.strand == Strand.NEG:
            nodes.reverse()
        return tuple(nodes)

    def _create_splice_graph(self):
        G = SGraph()
        for t in self.itertransfrags():
            nodes = [Exon(*n) for n in split_transfrag(t, self.node_bounds)]
            if self.strand == Strand.NEG:
                nodes.reverse()
            G.add_path(nodes, t.is_ref)
        # after graph is built, set node properties for start and stop nodes
        G.set_start_stop_nodes()
        return G

    def _mark_start_stop_nodes(self):
        G = self.G
        node_bounds = self.node_bounds
        strand = self.strand

        # mark start nodes
        changepts = set(self.start_sites)
        if self.guided_ends:
            changepts.update(self.ref_start_sites)
        for pos in changepts:
            if strand != Strand.NEG:
                bisect_func = bisect.bisect_right
            else:
                bisect_func = bisect.bisect_left
            i = bisect_func(node_bounds, pos)
            n = (node_bounds[i-1], node_bounds[i])
            if G.has_node(n):
                n_id = G.get_node_id(n)
                G.is_start[n_id] = True

        # mark stop nodes
        changepts = set(self.stop_sites)
        if self.guided_ends:
            changepts.update(self.ref_stop_sites)
        for pos in changepts:
            if strand == Strand.NEG:
                bisect_func = bisect.bisect_right
            else:
                bisect_func = bisect.bisect_left
            i = bisect_func(node_bounds, pos)
            n = (node_bounds[i-1], node_bounds[i])
            if G.has_node(n):
                n_id = G.get_node_id(n)
                G.is_stop[n_id] = True

    def _add_to_interval_trees(self, t):
        istart = Interval(t.start, t.start+1, value=t)
        iend = Interval(t.end-1, t.end, value=t)
        self.tree_starts.insert_interval(istart)
        self.tree_ends.insert_interval(iend)

    def _trim_change_point(self, cp):
        # search for matches in change point interval
        num_trimmed = 0
        if cp.sign < 0:
            hits = self.tree_ends.find(cp.pos, cp.end)
            for hit in hits:
                t = hit.value
                # last exon start cannot overlap interval
                last_exon_start = t.exons[-1][0]
                if cp.pos <= last_exon_start <= cp.end:
                    continue
                # trim end left
                t.exons[-1] = Exon(last_exon_start, cp.pos)
                num_trimmed += 1
        else:
            hits = self.tree_starts.find(cp.start, cp.pos)
            for hit in hits:
                t = hit.value
                # first exon end cannot overlap interval
                first_exon_end = t.exons[0][1]
                if cp.start <= first_exon_end <= cp.pos:
                    continue
                # trim start right
                t.exons[0] = Exon(cp.pos, first_exon_end)
                num_trimmed += 1

    def detect_change_points(self, *args, **kwargs):
        '''
        *args, **kwargs: passed directly to 'run_changepoint'

        returns list of ChangePoint tuples
        '''
        genome_id_str = ('%s:%d-%d[%s]' % (self.chrom, self.start, self.end,
                         Strand.to_gtf(self.strand)))
        changepts = []
        for n in self.G.nodes_iter():
            expr_data = self.get_expr_data(*n)
            for cp in run_changepoint(expr_data, *args, **kwargs):
                # add offset from start of node to change point positions
                cp = cp._replace(pos=n.start + cp.pos,
                                 start=n.start + cp.start,
                                 end=n.start + cp.end)
                changepts.append(cp)
                logging.debug('%s changepoint node=(%s-%s) '
                              'pos=%d interval=(%d-%d) p=%.3f fc=%.3f' %
                              (genome_id_str, n.start,
                               n.end, cp.pos, cp.start, cp.end,
                               cp.pvalue, cp.foldchange))
        return changepts

    def apply_change_point(self, cp, trim=True):
        '''
        trim: if true, will trim transfrags around change points
        '''
        if ((self.strand != Strand.NEG and cp.sign < 0) or
            (self.strand == Strand.NEG and cp.sign > 0)):
            self.stop_sites.add(cp.pos)
        else:
            self.start_sites.add(cp.pos)
        if trim:
            self._trim_change_point(cp)

    def get_change_point_gtf(self, cp):
        graph_id = ('G_%s_%d_%d_%s' %
                    (self.chrom, self.start, self.end,
                     Strand.to_gtf(self.strand)))
        features = []
        f = GTF.Feature()
        f.seqid = self.chrom
        f.source = 'taco'
        f.feature = 'changept'
        f.start = cp.pos
        f.end = cp.pos + 1
        f.score = 0
        f.strand = Strand.to_gtf(self.strand)
        f.phase = '.'
        f.attrs = {'graph_id': graph_id,
                   'sign': str(cp.sign),
                   'pvalue': str(cp.pvalue),
                   'foldchange': str(cp.foldchange)}
        features.append(f)
        f = GTF.Feature()
        f.seqid = self.chrom
        f.source = 'taco'
        f.feature = 'changeinterval'
        f.start = cp.start
        f.end = cp.end
        f.score = 0
        f.strand = Strand.to_gtf(self.strand)
        f.phase = '.'
        f.attrs = {'graph_id': graph_id,
                   'sign': str(cp.sign),
                   'pvalue': str(cp.pvalue),
                   'foldchange': str(cp.foldchange)}
        features.append(f)
        return features

    def split(self):
        '''splits into weakly connected component subgraphs'''
        # get connected components of graph which represent independent genes
        # unconnected components are considered different genes
        num_components, components = self.G.weakly_connected_components()
        if num_components == 1:
            yield self
            return
        # map transfrags to components
        component_transfrag_map = collections.defaultdict(list)
        # assign transfrags to components
        for t in self.itertransfrags():
            for n in split_transfrag(t, self.node_bounds):
                component_id = components[self.G.get_node_id(n)]
                component_transfrag_map[component_id].append(t)
                break
        # create new graphs using the separate components
        for component_transfrags in component_transfrag_map.itervalues():
            yield SpliceGraph.create(component_transfrags,
                                     guided_ends=self.guided_ends,
                                     guided_assembly=self.guided_assembly)

    def get_node_interval(self, n_id):
        return self.G.get_node(n_id)

    def get_node_id(self, n):
        return self.G.get_node_id(n)

    def get_node_expr_data(self, n_id):
        start, end = self.G.get_node(n_id)
        return self.get_expr_data(start, end)

    def itertransfrags(self):
        if self.guided_assembly:
            return chain(self.transfrags, self.ref_transfrags)
        else:
            return iter(self.transfrags)

    def get_start_stop_nodes(self):
        start_nodes = set()
        stop_nodes = set()
        for i in self.G.node_ids_iter():
            if self.G.is_start[i]:
                start_nodes.add(i)
            if self.G.is_stop[i]:
                stop_nodes.add(i)
        return start_nodes, stop_nodes

    def get_expr_data(self, start=None, end=None):
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if ((start < self.start) or (end > self.end)):
            m = ('query %d-%d outside locus bounds %d-%d' %
                 (start, end, self.start, self.end))
            raise TacoError(m)
        astart = start - self.start
        aend = end - self.start
        return self.expr_data[astart:aend]

    def get_node_gtf(self):
        graph_id = ('G_%s_%d_%d_%s' %
                    (self.chrom, self.start, self.end,
                     Strand.to_gtf(self.strand)))
        # iterate through locus and return change point data
        for n_id in self.G.node_ids_iter():
            n = self.get_node_interval(n_id)
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
            f.attrs = {'graph_id': graph_id,
                       'expr_min': str(expr_data.min()),
                       'expr_max': str(expr_data.max()),
                       'expr_mean': str(expr_data.mean()),
                       'ref_starts': ','.join(map(str, ref_starts)),
                       'ref_stops': ','.join(map(str, ref_stops))}
            yield f
