from array import array
from bisect import bisect_right
cimport cpython.array

# local imports
from sais cimport sais_int, sais_int_bwt

DEF TERMINATOR = 0
DEF SPACER = 1

cdef int suffix_cmp(tuple p, int[:] t, int[:] sa, int c):
    cdef int i = 0
    cdef int n = len(sa)
    cdef int np = len(p)
    cdef int retval = 0
    while i < len(p) and (sa[c]+i) < n:
        if p[i] > t[sa[c]+i]:
            retval = 1
            break
        elif p[i] < t[sa[c]+i]:
            retval = -1
            break
        i += 1
    return retval


cdef tuple suffix_range(tuple p, int[:] t, int[:] sa):
    cdef int n, l, r, s, c, retval

    n = len(sa)
    l = 0
    r = n
    while l < r:
        c = (l+r) // 2
        retval = suffix_cmp(p, t, sa, c)
        if retval > 0:
            l = c + 1
        else:
            r = c
    s = l
    r = n
    while l < r:
        c = (l+r) // 2
        retval = suffix_cmp(p, t, sa, c)
        if retval < 0:
            r = c
        else:
            l = c + 1
    return (s, r)


cdef class SuffixArrayIndex(object):

    cdef int[:] t
    cdef int[:] sa
    cdef list starts
    cdef list lengths
    cdef int alphabet_size

    def __init__(self, seqs):
        cdef object seq
        cdef list t, starts, lengths
        cdef int i, alphabet_size

        t = []
        starts = []
        lengths = []
        i = 0
        alphabet_size = 0
        for seq in seqs:
            t.extend(seq)
            t.append(SPACER)
            starts.append(i)
            lengths.append(len(seq))
            i += len(seq) + 1
        t.append(TERMINATOR)

        self.alphabet_size = max(t) + 1
        self.starts = starts
        self.lengths = lengths
        self.t = array('i', t)
        self.sa = array('i', t)
        if alphabet_size == 0:
            alphabet_size = max(t) + 1
        self.alphabet_size = alphabet_size
        suffix_array(self.t, self.sa, len(t), self.alphabet_size)

    def search(self, tuple p):
        cdef int i, l, r, si, istart, start, end

        l, r = suffix_range(p, self.t, self.sa)
        for i in xrange(l, r):
            si = self.sa[i]
            istart = bisect_right(self.starts, si) - 1
            start = self.starts[istart]
            end = start + self.lengths[istart]
            yield tuple(self.t[start:end])


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
