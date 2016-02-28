'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from cpython cimport array
import array

# local imports
from bsearch cimport bsearch_left, bsearch_right

def bisect_left(array.array a, int x, int lo=0, int hi=-1):
    if hi == -1:
        hi = len(a)
    return bsearch_left(a.data.as_ints, x, lo, hi)

def bisect_right(array.array a, int x, int lo=0, int hi=-1):
    if hi == -1:
        hi = len(a)
    return bsearch_right(a.data.as_ints, x, lo, hi)
