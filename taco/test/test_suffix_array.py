from taco.lib.suffixarray import SuffixArrayIndex
from taco.lib.csuffixarray import SuffixArrayIndex as CSuffixArrayIndex


def test_suffix_array():
    for IndexClass in (SuffixArrayIndex, CSuffixArrayIndex):
        kmers = [(1, 2), (1, 3), (1, 4), (1, 5)]
        sai = IndexClass(kmers)
        occ1 = list(sai.search((1,)))
        assert len(set(occ1).symmetric_difference(set(kmers))) == 0
        occ2 = list(sai.search((2,)))
        assert len(occ2) == 1 and (1, 2) in set(occ2)
        occ13 = list(sai.search((1, 3)))
        assert len(occ13) == 1 and (1, 3) in set(occ13)
        err1 = list(sai.search((-123,)))
        assert len(err1) == 0
        err1 = list(sai.search((9, 234)))
        assert len(err1) == 0


def test_different_length_kmers():
    for IndexClass in (SuffixArrayIndex, CSuffixArrayIndex):
        kmers = [(1,), (1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 2, 3, 4, 5),
                 (1, 3, 5), (1, 4, 5), (1, 2, 5)]
        sai = IndexClass(kmers)
        occ = set(sai.search((1,)))
        assert len(occ.symmetric_difference(set(kmers))) == 0
        occ = set(sai.search((2,)))
        assert len(occ) == 5 and (1, 2) in occ
        occ = set(sai.search((1, 2, 5)))
        assert len(occ) == 1 and (1, 2, 5) in occ
        occ = set(sai.search((2, 3, 4)))
        assert len(occ) == 2 and (1, 2, 3, 4, 5) in occ
