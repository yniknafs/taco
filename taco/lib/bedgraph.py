'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import collections
import numpy as np

from taco.lib.base import FLOAT_DTYPE

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

FMT_STRING = '%s\t%d\t%d\t%f\n'


def array_to_bedgraph(a, ref, start, fileh):
    if a.shape[0] == 0:
        return
    i = 0
    val = a[i]
    for j in xrange(1, a.shape[0]):
        newval = a[j]
        if val != newval:
            if val != 0:
                fileh.write(FMT_STRING % (ref, start + i, start + j, val))
            i = j
            val = newval
    if val != 0:
        fileh.write(FMT_STRING % (ref, start + i, start + a.shape[0], val))


def bedgraph_to_array(fileh):
    intervals = collections.defaultdict(lambda: [])
    maxend = collections.defaultdict(lambda: 0)
    for line in fileh:
        if line.startswith('#'):
            continue
        if not line:
            continue
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        ref, start, end, cov = (fields[0], int(fields[1]), int(fields[2]),
                                float(fields[3]))
        intervals[ref].append((start, end, cov))
        maxend[ref] = max(maxend[ref], end)
    covarrays = {}
    for ref, values in intervals.iteritems():
        covarrays[ref] = np.zeros(maxend[ref], dtype=FLOAT_DTYPE)
        for start, end, cov in values:
            covarrays[ref][start:end] += cov
    return covarrays
