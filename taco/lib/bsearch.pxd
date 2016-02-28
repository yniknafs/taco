# file: bsearch.pxd
cdef extern from "bsearch.h":
    int bsearch_right(int *a, int x, size_t lo, size_t hi);
    int bsearch_left(int *a, int x, size_t lo, size_t hi);
