'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import collections
import bisect
import numpy as np
import networkx as nx
import h5py

from base import Strand, TacoError, Exon
from dtypes import FLOAT_DTYPE, H5_CHUNKSIZE
from gtf import GTF
from cBedGraph import array_to_bedgraph
from cChangePoint import find_change_points

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def _find_boundaries(transfrags):
    chrom = None
    start = None
    end = None
    strands = set()
    for t in transfrags:
        if chrom is None:
            chrom = t.chrom
        elif chrom != t.chrom:
            raise TacoError('chrom mismatch')
        strands.add(t.strand)
        if start is None:
            start = t.start
        else:
            start = min(t.start, start)
        if end is None:
            end = t.end
        else:
            end = max(t.end, end)
    return chrom, start, end, strands


def _copy1d(a, astart, aend, bstart, bend):
    '''
    a: src array data for interval (astart, aend)
    b: dest array for interval (bstart, bend)
    '''
    asize = a.shape[0]
    bsize = bend - bstart
    assert asize == (aend - astart)
    if (bstart < astart) or (bend > aend):
        a1 = max(0, bstart - astart)
        a2 = min(asize, bend - astart)
        b1 = max(0, astart - bstart)
        b2 = min(bsize, aend - bstart)
        b = np.zeros((bend - bstart), dtype=FLOAT_DTYPE)
        b[b1:b2] = a[a1:a2]
        return b
    else:
        a1 = bstart - astart
        a2 = bend - astart
        return a[a1:a2]


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


class StrandedLocus(object):
    def __init__(self):
        self.guided_ends = False
        self.guided_assembly = False
        self.transfrags = []
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.expr_data = None
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        self.start_sites = set()
        self.stop_sites = set()

    @staticmethod
    def create(transfrags,
               guided_ends=False,
               guided_assembly=False):
        self = StrandedLocus()
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
        return self

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

    def create_splice_graph(self):
        '''
        returns networkx DiGraph object
        '''
        def add_node(G, n):
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=None)

        node_bounds = self._find_node_boundaries()
        G = nx.DiGraph()
        for t in self.transfrags:
            if t.is_ref and not self.guided_assembly:
                continue
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [Exon(a, b) for (a, b) in
                     split_transfrag(t, node_bounds)]
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
        G.graph['transfrags'] = self.transfrags
        G.graph['chrom'] = self.chrom
        G.graph['start'] = self.start
        G.graph['end'] = self.end
        G.graph['strand'] = self.strand
        return G

    def get_node_gtf(self):
        locus_id = ('L_%s_%d_%d_%s' %
                    (self.chrom, self.start, self.end,
                     Strand.to_gtf(self.strand)))
        # iterate through locus and return change point data
        G = self.create_splice_graph()
        for n in sorted(G):
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

    def split_splice_graph(self, G):
        # get connected components of graph which represent independent genes
        # unconnected components are considered different genes
        Gsubs = list(nx.weakly_connected_component_subgraphs(G))

        if len(Gsubs) == 1:
            Gsub = Gsubs[0]
            # set graph attributes
            Gsub.graph['transfrags'] = self.transfrags
            Gsub.graph['chrom'] = self.chrom
            Gsub.graph['start'] = self.start
            Gsub.graph['end'] = self.end
            Gsub.graph['strand'] = self.strand
            return Gsubs
        else:
            # map nodes to components
            node_subgraph_map = {}
            transfrag_subgraph_map = collections.defaultdict(list)
            for i, Gsub in enumerate(Gsubs):
                for n in Gsub:
                    node_subgraph_map[n] = i
            # assign transfrags to components
            for t in self.transfrags:
                if t.is_ref and not self.guided_assembly:
                    continue
                for n in split_transfrag(t, boundaries):
                    n = Exon(*n)
                    subgraph_id = node_subgraph_map[n]
                    transfrag_subgraph_map[subgraph_id].append(t)
                    break
            # set graph attributes
            for i, G in enumerate(Gsubs):
                # set graph attributes
                G.graph['transfrags'] = transfrag_subgraph_map[i]
                G.graph['chrom'] = self.chrom
                G.graph['start'] = self.start
                G.graph['end'] = self.end

                G.graph['strand'] = self.strand
            return Gsubs


class Locus(object):
    def __init__(self):
        self.guided_strand = False
        self.guided_ends = False
        self.guided_assembly = False
        self.chrom = None
        self.start = None
        self.end = None
        self.strands = set()
        self.expr_data = None
        self.strand_data = None
        self.strand_transfrags = [[], [], []]

    @staticmethod
    def create(transfrags,
               guided_strand=False,
               guided_ends=False,
               guided_assembly=False):
        self = Locus()
        self.guided_strand = guided_strand
        self.guided_ends = guided_ends
        self.guided_assembly = guided_assembly
        # find locus boundaries
        chrom, start, end, strands = _find_boundaries(transfrags)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strands.update(strands)
        # setup expression arrays
        self.expr_data = np.zeros((3, (end - start)), dtype=FLOAT_DTYPE)
        self.strand_data = np.zeros((3, (end - start)), dtype=np.bool)
        # add transcripts to locus
        for t in transfrags:
            self._add_transfrag(t)
        return self

    def _add_transfrag(self, t):
        if t.chrom != self.chrom:
            raise TacoError('chrom mismatch')
        if t.start < self.start:
            raise TacoError('transfrag start < locus start')
        if t.end > self.end:
            raise TacoError('transfrag end < locus end')
        self.strands.add(t.strand)
        self.strand_transfrags[t.strand].append(t)
        for exon in t.exons:
            astart = exon.start - self.start
            aend = exon.end - self.start
            if t.is_ref and self.guided_strand:
                self.strand_data[t.strand, astart:aend] = True
            else:
                self.expr_data[t.strand, astart:aend] += t.expr
                self.strand_data[t.strand, astart:aend] = True

    def create_stranded_loci(self):
        for s in self.strands:
            yield StrandedLocus.create(self.strand_transfrags[s],
                                       guided_ends=self.guided_ends,
                                       guided_assembly=self.guided_assembly)

    def _check_strand_ambiguous(self, t):
        '''
        Checks list of nodes for strandedness. If strand is unambiguous,
        then return pos or neg strand. If ambiguous, return unstranded.
        '''
        strands = [False, False]
        for exon in t.exons:
            start = exon.start - self.start
            end = exon.end - self.start
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
        unstranded_transfrags = self.strand_transfrags[Strand.NA]
        while(len(unstranded_transfrags) > 0):
            resolved = []
            unresolved = []
            for t in unstranded_transfrags:
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
            unstranded_transfrags = unresolved
            iterations += 1

        if num_resolved > 0:
            # clear unstranded locus
            self.strand_data[Strand.NA, :] = False
            self.expr_data[Strand.NA, :] = 0.0
            self.strands.remove(Strand.NA)
            # re-add unstranded transfrags
            for t in unstranded_transfrags:
                self._add_transfrag(t)
            self.strand_transfrags[Strand.NA] = unstranded_transfrags

        if iterations > 0:
            logging.debug('predict_unknown_strands: %d iterations' %
                          iterations)
        return num_resolved

    def get_expr_data(self, start, end, strand):
        if ((start < self.start) or (end > self.end)):
            m = ('query %d-%d outside locus bounds %d-%d' %
                 (start, end, self.start, self.end))
            raise TacoError(m)
        astart = start - self.start
        aend = end - self.start
        return self.expr_data[strand, astart:aend]

    def get_transfrags(self, strand):
        return self.strand_transfrags[strand]

    @staticmethod
    def open_bedgraphs(file_prefix):
        bgfilehs = []
        for s in (Strand.POS, Strand.NEG, Strand.NA):
            filename = '%s.%s.bedgraph' % (file_prefix, Strand.NAMES[s])
            bgfilehs.append(open(filename, 'w'))
        return bgfilehs

    @staticmethod
    def close_bedgraphs(bgfilehs):
        for fileh in bgfilehs:
            fileh.close()

    def write_bedgraph(self, bgfilehs):
        for s in self.strands:
            array_to_bedgraph(a=self.expr_data[s, :],
                              ref=self.chrom,
                              start=self.start,
                              fileh=bgfilehs[s])

    @staticmethod
    def open_expression_hdf5(filename, chrom_sizes_file):
        # read chrom sizes
        chrom_sizes = {}
        with open(chrom_sizes_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom_sizes[fields[0]] = int(fields[1])
        # create h5py datasets
        h5f = h5py.File(filename, mode='a', libver='latest')
        for ref, length in chrom_sizes.iteritems():
            chunksize = (3, min(H5_CHUNKSIZE, length))
            h5f.require_dataset(ref, shape=(3, length),
                                dtype=FLOAT_DTYPE,
                                exact=True,
#                                    chunks=True,
                                chunks=chunksize,
                                compression='lzf',
                                shuffle=True)
        return h5f

    def write_expression_hdf5(self, h5f):
        '''writes locus expression to 'h5f' an h5py.File object'''
        dset = h5f[self.chrom]
        dset[:, self.start:self.end] = self.expr_data
