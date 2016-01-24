'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import cStringIO
import timeit
import numpy as np

from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.bedgraph import array_to_bedgraph, bedgraph_to_array
from taco.lib.cBedGraph import array_to_bedgraph as c_array_to_bedgraph


def write_and_read_array(a, ref='chr1', start=0):
    buf = cStringIO.StringIO()
    array_to_bedgraph(a, ref, start, buf)
    contents = buf.getvalue()
    a = bedgraph_to_array(cStringIO.StringIO(contents))
    return a.get(ref, None)


def c_write_and_read_array(a, ref='chr1', start=0):
    filename = "tmp.bedgraph"
    with open(filename, 'w') as fileh:
        c_array_to_bedgraph(a, ref, start, fileh)
    a = bedgraph_to_array(open(filename))
    os.remove(filename)
    return a.get(ref, None)


def test_array1():
    a = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5],
                 dtype=FLOAT_DTYPE)
    y = write_and_read_array(a)
    assert np.array_equal(a, y)
    y = c_write_and_read_array(a)
    assert np.array_equal(a, y)


def test_array2():
    return
    a = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0], dtype=FLOAT_DTYPE)
    y = write_and_read_array(a)
    assert np.array_equal(a[:-3], y)
    y = c_write_and_read_array(a)
    assert np.array_equal(a[:-3], y)


def test_array3():
    return
    a = np.ones(5, dtype=FLOAT_DTYPE) * 10
    y = write_and_read_array(a)
    assert np.array_equal(a, y)
    y = c_write_and_read_array(a)
    assert np.array_equal(a, y)


def test_empty():
    return
    a = np.zeros(0, dtype=FLOAT_DTYPE)
    y = write_and_read_array(a)
    assert y is None
    y = c_write_and_read_array(a)
    assert y is None


def test_zeros():
    return
    a = np.zeros(5, dtype=FLOAT_DTYPE)
    y = write_and_read_array(a)
    assert y is None
    y = c_write_and_read_array(a)
    assert y is None


def test_performance():

    def stmt1():
        a = np.array(np.random.random(100000), dtype=FLOAT_DTYPE)
        buf = cStringIO.StringIO()
        array_to_bedgraph(a, ref='chr1', start=0, fileh=buf)
        # filename = "tmp.bedgraph"
        # with open(filename, 'w') as fileh:
        #     array_to_bedgraph(a, ref='chr1', start=0, fileh=fileh)
        # os.remove(filename)

    def stmt2():
        a = np.array(np.random.random(100000), dtype=FLOAT_DTYPE)
        filename = "tmp.bedgraph"
        with open(filename, 'w') as fileh:
            c_array_to_bedgraph(a, ref='chr1', start=0, fileh=fileh)
        os.remove(filename)

    t1 = timeit.Timer(stmt1)
    t2 = timeit.Timer(stmt2)
    print t1.timeit(number=2)
    print t2.timeit(number=2)

    # import pstats, cProfile
    # cProfile.runctx("stmt2()", globals(), locals(), "Profile.prof")
    # s = pstats.Stats("Profile.prof")
    # s.strip_dirs().sort_stats("time").print_stats()
