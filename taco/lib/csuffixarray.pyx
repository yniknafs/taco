cimport cpython.array

# local imports
from sais cimport sais_int, sais_int_bwt


def suffix_array(int[:] T, int[:] SA, int n, int k):
    sais_int(&T[0], &SA[0], n, k)


def suffix_array_bwt(int[:] T, int[:] U, int[:] A, int n, int k):
    sais_int_bwt(&T[0], &U[0], &A[0], n, k)


def bwt_from_sa(int[:] t, int[:] sa, end_value=0):
    '''Given T, returns BWT(T) by way of the suffix array'''
    cdef list bw
    cdef int si, end_row

    bw = []
    end_row = 0
    for si in sa:
        if si == 0:
            end_row = len(bw)
            bw.append(end_value)
        else:
            bw.append(t[si-1])
    return (bw, end_row)
