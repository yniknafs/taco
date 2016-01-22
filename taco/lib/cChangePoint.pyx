import numpy as np
cimport numpy as np
cimport cython

# import array type definitions
from dtypes import FLOAT_DTYPE
from dtypes cimport FLOAT_DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def find_change_points(np.ndarray[FLOAT_DTYPE_t, ndim=1] a,
                       int start=0,
                       FLOAT_DTYPE_t threshold=0):
    '''
    a - numpy array of expression data (all values >= 0)
    start - genomic start of 'a'
    threshold - expression threshold to detect changes
    '''
    cdef list change_points = []
    cdef int alen = a.shape[0]
    cdef int i

    if alen == 0:
        return []
    if alen == 1:
        return []
    for i in xrange(1, alen):
        if (a[i-1] > threshold) and (a[i] <= threshold):
            change_points.append(start + i)
        elif (a[i-1] <= threshold) and (a[i] > threshold):
            change_points.append(start + i)
    return change_points
