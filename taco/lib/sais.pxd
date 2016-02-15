# file: sais.pxd
cdef extern from "sais.h":
    int sais_int(const int *T, int *SA, int n, int k);
    int sais_int_bwt(const int *T, int *U, int *A, int n, int k);
