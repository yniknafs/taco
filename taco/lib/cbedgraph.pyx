import numpy as np
cimport numpy as np
cimport cython

from libc.stdio cimport FILE, fopen, fclose, fprintf, stdout

# define array types
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# for print statement
FMT_STRING = '%s\t%d\t%d\t%f\n'
DEF C_FMT_STRING = "%s\t%d\t%d\t%f\n"

@cython.boundscheck(False)
def array_to_bedgraph(np.ndarray[DTYPE_t, ndim=1] a,
                      char * ref, int start,
                      char * filename):
    assert a.dtype == DTYPE

    cdef int i, j, alen
    cdef DTYPE_t val, newval
    cdef FILE *fp

    alen = a.shape[0]
    if alen == 0:
        return

    # open output file
    fp = fopen(filename, 'w')

    i = 0
    val = a[i]
    for j in xrange(1, alen):
        newval = a[j]
        if val != newval:
            if val != 0:
                fprintf(fp, C_FMT_STRING, ref, start + i, start + j, val)
            i = j
            val = newval
    if val != 0:
        fprintf(fp, C_FMT_STRING, ref, start + i, start + alen, val)

    fclose(fp)
