'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import collections
import numpy as np
import h5py

from base import Strand, TacoError
from dtypes import FLOAT_DTYPE, H5_CHUNKSIZE
from cbedgraph import array_to_bedgraph
from splice_graph import SpliceGraph

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

    def create_splice_graphs(self, split=True):
        for s in self.strands:
            sgraph = SpliceGraph.create(self.strand_transfrags[s],
                                        guided_ends=self.guided_ends,
                                        guided_assembly=self.guided_assembly)
            if not split:
                yield sgraph
            else:
                for ssubgraph in sgraph.split():
                    yield ssubgraph

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
    def get_bedgraph_file_names(file_prefix):
        for s in (Strand.POS, Strand.NEG, Strand.NA):
            filename = '%s.%s.bedgraph' % (file_prefix, Strand.NAMES[s])
            yield (s, filename)

    @staticmethod
    def open_bedgraphs(file_prefix):
        bgfilehs = []
        for s, filename in Locus.get_bedgraph_file_names(file_prefix):
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
    def open_splice_bed(filename):
        fh = open(filename, 'w')
        track_line = ('track name=junctions description="Splice Junctions" '
                      'graphType=junctions')
        print >>fh, track_line
        return fh

    def write_splice_bed(self, fh):
        intron_dict = collections.defaultdict(float)
        for strand in Strand.POS, Strand.NEG:
            for t in self.strand_transfrags[strand]:
                if t.is_ref:
                    continue
                for start, end in t.iterintrons():
                    intron_dict[(start, end, strand)] += t.expr
        for intron, expr in intron_dict.iteritems():
            start, end, strand = intron
            fields = [self.chrom, str(start - 1), str(end + 1), 'JUNC',
                      str(expr), Strand.to_gtf(strand),
                      str(start - 1), str(end + 1), '255,0,0',
                      '2', '1,1', '0,%d' % (end + 1 - start)]
            print >>fh, '\t'.join(fields)
