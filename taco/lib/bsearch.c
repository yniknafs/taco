#include <stdio.h>

/*
 * bsearch_right - binary search an array of integers
 * a: pointer to int array
 * n: size of array
 * x: key to search for
 */
int bsearch_right(int *a, int x, size_t lo, size_t hi)
{
    size_t mid;
    int res;

    while (lo < hi) {
        mid = ((size_t)lo + hi) / 2;
        res = x - a[mid];
        if (res < 0)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

/*
 * bsearch_left - binary search an array of integers
 * a: pointer to int array
 * n: size of array
 * x: key to search for
 */
int bsearch_left(int *a, int x, size_t lo, size_t hi)
{
    size_t mid;
    int res;

    while (lo < hi) {
        mid = ((size_t)lo + hi) / 2;
        res = x - a[mid];
        if (res <= 0)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

int main ()
{
    int values[] = { 5, 20, 29, 32, 32, 32, 63 };
    int n = 7;
    int key = 32;
    int index = 0;

    index = bsearch_left(values, key, 0, n);
    printf("Index=%d\n", index);
    index = bsearch_right(values, key, 0, n);
    printf("Index=%d\n", index);

    return 0;
}
