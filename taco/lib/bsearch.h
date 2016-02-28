#ifndef _BSEARCH_H
#define _BSEARCH_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int bsearch_right(int *a, int x, size_t lo, size_t hi);
int bsearch_left(int *a, int x, size_t lo, size_t hi);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _BSEARCH_H */
