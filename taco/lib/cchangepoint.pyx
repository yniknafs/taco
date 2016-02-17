import numpy as np
cimport numpy as np
cimport cython

# import array type definitions
from dtypes import FLOAT_DTYPE
from dtypes cimport FLOAT_DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def find_threshold_points(np.ndarray[FLOAT_DTYPE_t, ndim=1] a,
                          int start=0, FLOAT_DTYPE_t threshold=0):
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

@cython.boundscheck(False)
@cython.wraparound(False)
def mse(np.ndarray[FLOAT_DTYPE_t, ndim=1] x):
    cdef int xsize = x.shape[0]
    if xsize < 1:
        return 0.0, -1
    cdef FLOAT_DTYPE_t x1sum, x2sum, x1mean, x2mean
    cdef FLOAT_DTYPE_t mse, mse_min
    cdef int mse_i = -1
    cdef int i, j, prev

    cdef list changepts = []
    for i in xrange(1, xsize):
        if x[i] != x[i-1]:
            changepts.append(i)
    if len(changepts) == 0:
        return 0.0, -1

    x1sum = 0
    x2sum = x.sum()
    mse_i = -1
    prev = 0
    for i in changepts:
        # compute new means
        for j in xrange(prev, i):
            x1sum += x[j]
            x2sum -= x[j]
        x1mean = x1sum / i
        x2mean = x2sum / (xsize - i)

        # within cluster sum of squares
        mse = 0.0
        for j in xrange(0, i):
            mse += (x[j] - x1mean) ** 2
        for j in xrange(i, xsize):
            mse += (x[j] - x2mean) ** 2

        if (mse_i == -1) or (mse < mse_min):
            mse_min = mse
            mse_i = i
        prev = i
    return mse_min, mse_i
