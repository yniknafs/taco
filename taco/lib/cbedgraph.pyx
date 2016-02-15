'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np
cimport numpy as np
cimport cython

# import array type definitions
from dtypes import FLOAT_DTYPE
from dtypes cimport FLOAT_DTYPE_t

# declare the Python macro to access files:
cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    void fprintf(FILE* f, char *s, ...)

# declare builtin file object
cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

# for print statement
FMT_STRING = '%s\t%d\t%d\t%f\n'
DEF C_FMT_STRING = "%s\t%d\t%d\t%f\n"

# Now declare the C function that requires a file:
cdef void c_write(FILE* fp, char* ref, int start, int end, FLOAT_DTYPE_t val):
    fprintf(fp, C_FMT_STRING, ref, start, end, val)


@cython.boundscheck(False)
def array_to_bedgraph(np.ndarray[FLOAT_DTYPE_t, ndim=1] a,
                      char * ref, int start,
                      file fileh):
    assert a.dtype == FLOAT_DTYPE

    cdef int i, j, alen
    cdef FLOAT_DTYPE_t val, newval
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
