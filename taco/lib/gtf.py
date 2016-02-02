'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer
'''
import os
import collections
import subprocess

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2015"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def sort_gtf(filename, output_file, tmp_dir=None):
    args = ["sort"]
    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", "-k4,4n", "-k3,3r", filename])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)


class GTFError(Exception):
    pass


class GTF:
    POS_STRAND = '+'
    NEG_STRAND = '-'
    NO_STRAND = '.'
    EMPTY_FIELD = '.'
    ATTR_SEP = '; '

    class Attr:
        GENE_ID = 'gene_id'
        TRANSCRIPT_ID = 'transcript_id'
        SAMPLE_ID = 'sample_id'
        REF = 'ref'
        EXPR = 'expr'

    class Feature:
        '''
        GTF Specification (fields are tab-separated)
        1. seqname - sequence name (chromosome)
        2. source - program that generated this feature.
        3. feature - type of feature ("transcript", "exon")
        4. start - start pos of feature (1-based)
        5. end - end pos of feature (inclusive)
        6. score - number between 0 and 1000
        7. strand - '+', '-', or '.' for unspecified
        8. phase - for coding exons, frame should be a number between 0-2
        that represents the reading frame of the first base. If the feature
        is not a coding exon, the value should be '.'.
        9. attributes - optional attributes in the format:
            key1 "value1"; key2 "value2";
        '''
        __slots__ = ('seqid', 'source', 'feature', 'start', 'end',
                     'score', 'strand', 'phase', 'attrs')

        def __init__(self):
            pass

        def __str__(self):
            attr_str = ' '.join('%s "%s";' % (k, v)
                                for (k, v) in self.attrs.iteritems())
            fields = [self.seqid,
                      self.source,
                      self.feature,
                      # convert to 1-based intervals
                      str(self.start + 1),
                      str(self.end),
                      str(self.score),
                      self.strand,
                      self.phase,
                      attr_str]
            return '\t'.join(fields)

        @staticmethod
        def from_str(s):
            fields = s.strip().split('\t')
            f = GTF.Feature()
            f.seqid = fields[0]
            f.source = fields[1]
            f.feature = fields[2]
            # convert from 1-based (inclusive) to 0-based (exclusive) intervals
            f.start = int(fields[3])-1
            f.end = int(fields[4])
            f.score = fields[5]
            f.strand = fields[6]
            f.phase = fields[7]
            attrs = collections.OrderedDict()  # preserve readability
            if fields[8] != GTF.EMPTY_FIELD:
                attr_strings = fields[8].strip().split(';')
                for a in attr_strings:
                    if not a.strip():
                        continue
                    # key-value pair separated by whitespace
                    k, v = a.strip().split()
                    # remove quotes
                    v = v.strip('"')
                    attrs[k] = v
            f.attrs = attrs
            return f

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            yield GTF.Feature.from_str(line)

    @staticmethod
    def parse_loci(line_iter):
        '''
        Parses appropriately formatted GTF file into individual 'loci'. A
        locus is defined as a set of overlapping features.

        Assumes that GTF file has been sorted and formatted such that a
        single 'transcript' feature appears before its corresponding 'exon'
        features (this simplifies parsing).
        '''
        def window_overlap(a, b):
            if a[0] != b[0]:
                return False
            # return (a[1] <= b[2]) and (b[1] <= a[2])
            return (a[1] < b[2]) and (b[1] < a[2])

        def get_intervals(line_iter):
            for line in line_iter:
                if line.startswith("#"):
                    continue
                # read the essential part of the GTF line
                line = line.rstrip()
                fields = line.split('\t', 5)
                if len(fields) < 5:
                    continue
                seqid = fields[0]
                start = int(fields[3])-1
                end = int(fields[4])
                yield seqid, start, end, line

        try:
            interval_iter = get_intervals(line_iter)
            # initialize window
            seqid, start, end, line = interval_iter.next()
            window = [line]
            window_range = (seqid, start, end)
            # separate into loci
            for seqid, start, end, line in interval_iter:
                # check if next transcript is outside current window
                interval = (seqid, start, end)
                if not window_overlap(interval, window_range):
                    # yield current window
                    yield window_range, window
                    # reset window
                    window = [line]
                    window_range = (seqid, start, end)
                else:
                    # add transcript to window
                    window.append(line)
                    newstart = (start if start < window_range[1]
                                else window_range[1])
                    newend = (end if end > window_range[2]
                              else window_range[2])
                    window_range = (seqid, newstart, newend)
        except StopIteration:
            pass

        # yield last window
        if len(window) > 0:
            yield window_range, window
