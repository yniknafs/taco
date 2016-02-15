from array import array
from itertools import chain
from bisect import bisect_left, bisect_right

from taco.lib.csuffixarray import suffix_array
from taco.lib.fmindex import FmIndex


def test_same_length_kmers():
    kmers = [(2, 3, 4, 5, 6), (1, 2, 3, 4, 5), (3, 4, 5, 6, 7),
             (5, 6, 7, 8, 9), (4, 5, 6, 7, 8), (6, 7, 8, 9, 10)]
    t = list(chain(*kmers))
    t.append(0)
    t = array('i', t)
    fm = FmIndex(t, alphabet_size=256)
    result = set()
    for i in fm.occurrences((3, 4, 5)):
        start = 5 * (i / 5)
        end = start + 5
        result.add(tuple(t[start:end]))
    assert len(result) == 3
    assert (2, 3, 4, 5, 6) in result
    assert (1, 2, 3, 4, 5) in result
    assert (3, 4, 5, 6, 7) in result


def test_different_length_kmers():
    kmers = [(1,), (2, 3), (3, 4, 5), (2, 3, 4, 5, 6), (3, 4, 5, 6, 7),
             (4, 5, 6), (5, 6, 7, 8, 9), (6, 7, 8, 9, 10), (1, 3, 4, 5, 10),
             (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]
    t = []
    starts = []
    lengths = []
    i = 0
    for kmer in kmers:
        t.extend(kmer)
        starts.append(i)
        lengths.append(len(kmer))
        i += len(kmer)
    t.append(0)
    t = array('i', t)
    fm = FmIndex(t, alphabet_size=256)

    print 't', t
    print 'len t', len(t)

    for i in fm.occurrences((3, 4, 5)):
        print 'i', i
        istart = bisect_right(starts, i) - 1
        start = starts[istart]
        length = lengths[istart]
        print 'yo result', t[start:start+length]

    print '----'
    for i in fm.occurrences((1,)):
        istart = bisect_right(starts, i) - 1
        start = starts[istart]
        length = lengths[istart]
        print 'result', t[start:start+length]

    return

    result = set()
    for i in fm.occurrences((3, 4, 5)):
        start = 5 * (i / 5)
        end = start + 5
        result.add(tuple(t[start:end]))
    assert len(result) == 3
    assert (2, 3, 4, 5, 6) in result
    assert (1, 2, 3, 4, 5) in result
    assert (3, 4, 5, 6, 7) in result
