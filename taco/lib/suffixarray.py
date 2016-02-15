'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from array import array
from bisect import bisect_right

from csuffixarray import suffix_array

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def _cmp_suffix(p, t, sa, c, n):
    pgt = False
    plt = False
    i = 0
    while i < len(p) and (sa[c]+i) < n:
        if p[i] > t[sa[c]+i]:
            pgt = True
            break
        elif p[i] < t[sa[c]+i]:
            plt = True
            break
        i += 1
    return pgt, plt


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
        n = len(t)
        l, r = 0, n
        while l < r:
            c = (l+r) // 2
            pgt, plt = _cmp_suffix(p, t, sa, c, n)
            if pgt:
                l = c + 1
            else:
                r = c
        s = l
        r = n
        while l < r:
            c = (l+r) // 2
            pgt, plt = _cmp_suffix(p, t, sa, c, n)
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
            yield tuple(t[start:end])
