'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import numpy as np

from base import Strand, TacoError, FLOAT_DTYPE
from cbedgraph import array_to_bedgraph

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def find_genomic_boundaries(transfrags):
    # determine genomic bounds of all transfrags
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
        start = t.start if start is None else min(t.start, start)
        end = t.end if end is None else max(t.end, end)
    return chrom, start, end, strands


class Locus(object):

    def __init__(self):
        self.chrom = None
        self.start = None
        self.end = None
        self.strands = None
        self.expr_data = None
        self.strand_data = None
        self.strand_transfrags = [[], [], []]
        self.guided_strand = False

    @staticmethod
    def create(transfrags, guided_strand=False):
        self = Locus()
        self.guided_strand = guided_strand
        # define locus boundaries
        chrom, start, end, strands = find_genomic_boundaries(transfrags)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strands = strands
        # initialize expression data arrays
        length = self.end - self.start
        self.expr_data = np.zeros((3, length), dtype=FLOAT_DTYPE)
        self.strand_data = np.zeros((3, length), dtype=np.bool)
        # add transcripts to locus
        for t in transfrags:
            self.strand_transfrags[t.strand].append(t)
            self._add_transfrag(t)
        return self

    def _add_transfrag(self, t):
        for exon in t.exons:
            start = exon.start - self.start
            end = exon.end - self.start
            self.expr_data[t.strand, start:end] += t.expr
            if ((not t.is_ref) or (t.is_ref and self.guided_strand)):
                self.strand_data[t.strand, start:end] = True

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
        if ((start < self.start) or (end > self.end)):
            raise TacoError('start/end are out of bounds')
        start = start - self.start
        end = end - self.start
        return self.expr_data[strand, start:end]

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
        for strand in (Strand.POS, Strand.NEG, Strand.NA):
            array_to_bedgraph(a=self.expr_data[strand],
                              ref=self.chrom,
                              start=self.start,
                              fileh=bgfilehs[strand])
