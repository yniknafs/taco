'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from cpython cimport array
cimport cbisect

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


cdef list split_intervals(int *starts, int *ends, size_t n, int *boundaries):
    cdef list nodes
    cdef int i, start_ind, end_ind
    nodes = []
    for i in xrange(n):
        start_ind = cbisect.bisect_right(boundaries, starts[i], lo=end_ind)
        end_ind = cbisect.bisect_left(boundaries, ends[i], lo=start_ind)
        if start_ind == end_ind:
            nodes.append((boundaries[start_ind] - 1, boundaries[start_ind]))
        else:
            for j in xrange(start_ind-1, end_ind):
                nodes.append((boundaries[j], boundaries[j+1]))
    return nodes


def split_transfrag(list exons, array.array boundaries):
    '''
    exons must be sorted in increasing order

    output: (generator) tuples (a,b) reflecting nodes
    '''
    cdef int i, start_ind, end_ind
