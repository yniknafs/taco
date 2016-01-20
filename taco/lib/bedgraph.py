'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def array_to_bedgraph(a, omit_zeros=True):
    if a.shape[0] == 0:
        return
    i = 0
    val = a[i]
    for j in xrange(1, a.shape[0]):
        newval = a[j]
        if val != newval:
            if not (omit_zeros and (val == 0)):
                yield i, j, val
            i = j
            val = newval
    if i != a.shape[0]:
        if not (omit_zeros and (val == 0)):
            yield i, a.shape[0], val
