'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import bisect
import numpy as np
import networkx as nx
import h5py

from base import Strand, TacoError, Exon
from dtypes import FLOAT_DTYPE, H5_CHUNKSIZE
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
            print boundaries, start_ind, end_ind, exon.start, exon.end, 'isref', t.is_ref, 'SDFLKJSDF'
            yield boundaries[start_ind-1], boundaries[start_ind]
        else:
            # all the splice sites in between the exon borders must overlap
            for j in xrange(start_ind-1, end_ind):
                yield boundaries[j], boundaries[j+1]


class StrandedLocus(object):
    def __init__(self, chrom, start, end, strand,
                 guided_strand=False,
                 guided_ends=False,
                 guided_assembly=False):
        self.guided_strand = guided_strand
        self.guided_ends = guided_ends
        self.guided_assembly = guided_assembly
        self.transfrags = []
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.expr_data = np.zeros((end - start), dtype=FLOAT_DTYPE)
        self.strand_data = np.zeros((end - start), dtype=np.bool)
        self.ref_start_sites = set()
        self.ref_stop_sites = set()
        self.start_sites = set()
        self.stop_sites = set()
        self.node_bounds = None

    def add_transfrag(self, t):
        if t.chrom != self.chrom:
            raise TacoError('chrom mismatch')
        if t.strand != self.strand:
            raise TacoError('strand mismatch')
        if t.start < self.start:
            raise TacoError('transfrag start < locus start')
        if t.end > self.end:
            raise TacoError('transfrag end < locus end')
        self.transfrags.append(t)
        for exon in t.exons:
            astart = exon.start - self.start
            aend = exon.end - self.start
            if t.is_ref:
                self.ref_start_sites.add(t.txstart)
                self.ref_stop_sites.add(t.txstop)
                if self.guided_strand:
                    self.strand_data[astart:aend] = True
            else:
                self.expr_data[astart:aend] += t.expr
                self.strand_data[astart:aend] = True

    def get_expr_data(self, start, end):
        return _copy1d(self.expr_data, self.start, self.end,
                       start, end)

    def write_bedgraph(self, fileh):
        array_to_bedgraph(a=self.expr_data,
                          ref=self.chrom,
                          start=self.start,
                          fileh=fileh)

    def write_expression_hdf5(self, h5f):
        '''writes locus expression to 'h5f' an h5py.File object'''
        grp = '%s/%s' % (self.chrom, Strand.NAMES[self.strand])
        dset = h5f[grp]
        dset[self.start:self.end] = self.expr_data

    def _find_node_boundaries(self):
        node_bounds = set()
        # nodes bounded by regions where expression changes to/from zero
        zero_points = find_change_points(self.expr_data, start=self.start)
        node_bounds.update(zero_points)
        # nodes bounded by introns
        start = None
        end = None
        for t in self.transfrags:
            if t.is_ref:
                if self.guided_ends:
                    node_bounds.update((t.start, t.end))
                if self.guided_assembly:
                    node_bounds.update(t.itersplices())
            else:
                node_bounds.update(t.itersplices())
                if start is None:
                    start = t.start
                else:
                    start = min(t.start, start)
                if end is None:
                    end = t.end
                else:
                    end = max(t.end, end)
        node_bounds.add(start)
        node_bounds.add(end)
        print 'node bounds', sorted(node_bounds)
        return sorted(node_bounds)

    def create_splice_graph(self):
        def add_node(G, n):
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=None)

        boundaries = self._find_node_boundaries()
        G = nx.DiGraph()
        for t in self.transfrags:
            if t.is_ref and not self.guided_assembly:
                continue
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [Exon(a, b) for (a, b) in
                     split_transfrag(t, boundaries)]
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
        G.graph['boundaries'] = boundaries
        G.graph['chrom'] = self.chrom
        G.graph['start'] = boundaries[0]
        G.graph['end'] = boundaries[-1]
        G.graph['strand'] = self.strand
        G.graph['expr_data'] = self.get_expr_data(boundaries[0], boundaries[1])
        return G


class Locus(object):
    def __init__(self):
        self.chrom = None
        self.start = None
        self.end = None
        self.strands = None
        self.L = [None, None, None]

    @staticmethod
    def create(transfrags,
               guided_strand=False,
               guided_ends=False,
               guided_assembly=False):
        self = Locus()
        # find locus boundaries
        chrom, start, end, strands = _find_boundaries(transfrags)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strands = strands
        # create stranded loci
        for s in self.strands:
            self.L[s] = StrandedLocus(chrom=chrom,
                                      start=start,
                                      end=end,
                                      strand=s,
                                      guided_strand=guided_strand,
                                      guided_ends=guided_ends,
                                      guided_assembly=guided_assembly)
        # add transcripts to locus
        for t in transfrags:
            slocus = self.L[t.strand]
            slocus.add_transfrag(t)
        return self

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
                locus = self.L[s]
                if locus and locus.strand_data[start:end].any():
                    strands[s] = True
                if locus and locus.expr_data[start:end].sum() > 0:
                    strands[s] = True
        if strands[Strand.POS] and strands[Strand.NEG]:
            return Strand.NA
        elif strands[Strand.POS]:
            return Strand.POS
        elif strands[Strand.NEG]:
            return Strand.NEG
        return Strand.NA

    def impute_unknown_strands(self):
        if self.L[Strand.NA] is None:
            return
        unstranded_transfrags = self.L[Strand.NA].transfrags
        # iteratively predict strand of unstranded transfrags
        # stop when no new transfrag strands can be imputed
        iterations = 0
        num_resolved = 0
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
                self.L[new_strand].add_transfrag(t)
            unstranded_transfrags = unresolved
            iterations += 1

        if num_resolved > 0:
            # clear unstranded locus and re-add unresolved transcripts
            oldlocus = self.L[Strand.NA]
            locus = StrandedLocus(chrom=oldlocus.chrom,
                                  start=oldlocus.start,
                                  end=oldlocus.end,
                                  strand=oldlocus.strand,
                                  guided_strand=oldlocus.guided_strand,
                                  guided_ends=oldlocus.guided_ends,
                                  guided_assembly=oldlocus.guided_assembly)
            for t in unstranded_transfrags:
                locus.add_transfrag(t)
            self.L[Strand.NA] = locus

        if iterations > 0:
            logging.debug('predict_unknown_strands: %d iterations' %
                          iterations)
        return num_resolved

    def get_expr_data(self, start, end, strand):
        if ((start < self.start) or (end > self.end)):
            m = ('query %d-%d outside locus bounds %d-%d' %
                 (start, end, self.start, self.end))
            raise TacoError(m)
        locus = self.L[strand]
        if locus is None:
            return np.zeros(end - start, dtype=FLOAT_DTYPE)
        return locus.get_expr_data(start, end)

    def get_transfrags(self, strand):
        slocus = self.L[strand]
        if slocus is None:
            return []
        return slocus.transfrags

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
        for s in (Strand.POS, Strand.NEG, Strand.NA):
            slocus = self.L[s]
            if slocus is None:
                continue
            slocus.write_bedgraph(bgfilehs[s])

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
            for s in Strand.STRANDS:
                grp = '%s/%s' % (ref, Strand.NAMES[s])
                chunksize = min(H5_CHUNKSIZE, length)
                h5f.require_dataset(grp, shape=(length,),
                                    dtype=FLOAT_DTYPE,
                                    exact=True,
#                                    chunks=True,
                                    chunks=(chunksize,),
                                    compression='lzf',
                                    shuffle=True)
        return h5f

    def write_expression_hdf5(self, h5f):
        '''writes locus expression to 'h5f' an h5py.File object'''
        for slocus in self.L:
            if slocus is None:
                continue
            slocus.write_expression_hdf5(h5f)
