'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import random
from array import array
from bisect import bisect_left, bisect_right

import taco.lib.cbisect as cbisect


def test_bisect():
    a = array('i', [5, 20, 29, 32, 32, 32, 63])
    for x in xrange(0, max(a) + 1):
        assert cbisect.bisect_left(a, x) == bisect_left(a, x)
        assert cbisect.bisect_right(a, x) == bisect_right(a, x)
    a = array('i', sorted([random.randint(0, 1000) for x in xrange(10000)]))
    for x in xrange(0, max(a) + 1):
        assert cbisect.bisect_left(a, x) == bisect_left(a, x)
        assert cbisect.bisect_right(a, x) == bisect_right(a, x)
