from array import array
from itertools import chain
from bisect import bisect_left, bisect_right

from taco.lib.csuffixarray import suffix_array
from taco.lib.fmindex import FmIndex


def find3(L, S):
    # Suffix array
    start = 0
    end = len(SL)
    while start < end:
        mid = start + (end - start) // 2
        pa = SL_key_fn(W, SL[mid], 100)
        pb = SL_key_fn(S, 0, len(S))
        if pa < pb:
            start = mid + 1
        elif pb < pa:
            end = mid
        else:
            # A word may contain the same S multiple times
            R = set()
            while mid > 0 and W.startswith(S, SL[mid]):
                mid = mid - 1
            if not W.startswith(S, SL[mid]):
                mid = mid + 1
            while mid < len(SL) and W.startswith(S, SL[mid]):
                p = bisect.bisect_right(T, SL[mid]) - 1
                e = L[p]
                assert S in e
                R.add(p)
                mid = mid + 1
            return [(L[i]) for i in R]
    return []


def binarySearchSA(t, sa, p):
    assert t[-1] == 0 # t already has terminator
    assert len(t) == len(sa) # sa is the suffix array for t
    if len(t) == 1: return 1
    l, r = 0, len(sa) # invariant: sa[l] < p < sa[r]
    while True:
        c = (l + r) // 2
        # determine whether p < T[sa[c]:] by doing comparisons
        # starting from left-hand sides of p and T[sa[c]:]
        plt = True # assume p < T[sa[c]:] until proven otherwise
        i = 0
        while i < len(p) and sa[c]+i < len(t):
            if p[i] < t[sa[c]+i]:
                break # p < T[sa[c]:]
            elif p[i] > t[sa[c]+i]:
                plt = False
                break # p > T[sa[c]:]
            i += 1 # tied so far
        if plt:
            if c == l + 1: return c
            r = c
        else:
            if c == r - 1: return r
            l = c


def suffix_array_binary_search(t, sa, p):
    assert t[-1] == 0
    assert len(t) == len(sa)
    l, r = 0, len(sa)
    while l < r:
        c = (l+r) // 2
        pgt = True
        i = 0
        while i < len(p) and sa[c]+i < len(t):
            if p[i] > t[sa[c]+i]:
                break
            else:
                pgt = False
                break
            i += 1
        if pgt:
            l = c + 1
        else:
            r = c
    s = l
    r = len(sa)
    while l < r:
        c = (l+r) // 2
        plt = True
        i = 0
        while i < len(p) and sa[c]+i < len(t):
            if p[i] < t[sa[c]+i]:
                break
            else:
                plt = False
                break

            i += 1
        if plt:
            r = c
        else:
            l = c + 1
    return (s, r)


class SuffixArrayIndex(object):

    def __init__(self, seqs, alphabet_size=0):
        t = []
        starts = []
        lengths = []
        i = 0
        for seq in seqs:
            t.extend(seq)
            starts.append(i)
            lengths.append(len(seq))
            i += len(seq)
        t.append(0)
        self.starts = starts
        self.lengths = lengths
        self.t = array('i', t)
        self.sa = array('i', t)
        if alphabet_size == 0:
            alphabet_size = len(set(t))
        self.alphabet_size = max(alphabet_size, 256)
        suffix_array(self.t, self.sa, len(t), self.alphabet_size)

    def _range(self, p):
        t = self.t
        sa = self.sa
        assert t[-1] == 0
        assert len(t) == len(sa)
        l, r = 0, len(sa)
        while l < r:
            c = (l+r) // 2
            pgt = True
            i = 0
            while i < len(p) and sa[c]+i < len(t):
                if p[i] > t[sa[c]+i]:
                    break
                else:
                    pgt = False
                    break
                i += 1
            if pgt:
                l = c + 1
            else:
                r = c
        s = l
        r = len(sa)
        while l < r:
            c = (l+r) // 2
            plt = True
            i = 0
            while i < len(p) and sa[c]+i < len(t):
                if p[i] < t[sa[c]+i]:
                    break
                else:
                    plt = False
                    break

                i += 1
            if plt:
                r = c
            else:
                l = c + 1
        return (s, r)

    def search(self, p):
        starts = self.starts
        lengths = self.lengths
        t = self.t
        sa = self.sa
        l, r = self._range(p)
        for i in xrange(l, r):
            si = sa[i]
            istart = bisect_right(starts, si) - 1
            start = starts[istart]
            end = start + lengths[istart]
            yield t[start:end]


def test_suffix_array():
    pass


def test_suffix_array():
    kmers = [(2, 3, 4, 5, 6), (1, 2, 3, 4, 5), (3, 4, 5, 6, 7),
             (5, 6, 7, 8, 9), (4, 5, 6, 7, 8), (6, 7, 8, 9, 10),
             (1, 2, 3)]
    t = list(chain(*kmers))
    t.append(0)
    t = array('i', t)
    sa = array('i', t)
    suffix_array(t, sa, len(t), 256)


    sai = SuffixArrayIndex(kmers)
    print list(sai.search((1,)))

    return

    a, b = suffix_array_binary_search(t, sa, (1,))
    print 'occ', a, b
    for occ in xrange(a, b):
        print 'occ', occ, 'sa[occ]', sa[occ], 't[sa[occ]]', t[sa[occ]]

    print 'search', binarySearchSA(t, sa, (1,))
    print 'search2', suffix_array_binary_search(t, sa, (1,))
    print '----'
    print 't', t
    print 'sa', sa
    print '-----'
    return
    a, b = suffix_array_binary_search(t, sa, (1,))
    print 'hi'
    print a, b
    print 'bye'



def test_same_length_kmers():
    return
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
    return
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
