import numpy as np
cimport numpy as np
cimport cython

# declare the Python macro to access files:
cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    void fprintf(FILE* f, char *s, ...)


# declare builtin file object
cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

# define array types
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# for print statement
FMT_STRING = '%s\t%d\t%d\t%f\n'
DEF C_FMT_STRING = "%s\t%d\t%d\t%f\n"

# Now declare the C function that requires a file:
cdef void c_write(FILE* fp, char* ref, int start, int end, DTYPE_t val):
    fprintf(fp, C_FMT_STRING, ref, start, end, val)



@cython.boundscheck(False)
def array_to_bedgraph(np.ndarray[DTYPE_t, ndim=1] a,
                      char * ref, int start,
                      file fileh):
    assert a.dtype == DTYPE

    cdef int i, j, alen
    cdef DTYPE_t val, newval
    cdef FILE *fp = PyFile_AsFile(fileh)

    alen = a.shape[0]
    if alen == 0:
        return

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
