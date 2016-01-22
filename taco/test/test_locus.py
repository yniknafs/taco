'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import Strand, TacoError
from taco.lib.transfrag import Transfrag
from taco.lib.locus import Locus, _copy1d
from taco.test.base import read_gtf


def test_copy1d():
    a = np.array([1, 2, 3, 4, 5])
    astart, aend = 100, 105
    # perfect copy
    b = _copy1d(a, astart, aend, 100, 105)
    assert np.array_equal(b, a)
    # copy outside left bounds
    b = _copy1d(a, astart, aend, 90, 100)
    assert np.array_equal(b, np.zeros(10))
    # copy outside right bounds
    b = _copy1d(a, astart, aend, 105, 115)
    assert np.array_equal(b, np.zeros(10))
    # cross into left bound
    b = _copy1d(a, astart, aend, 98, 102)
    assert np.array_equal(b, [0, 0, 1, 2])
    # cross into right bound
    b = _copy1d(a, astart, aend, 103, 107)
    assert np.array_equal(b, [4, 5, 0, 0])
    # outside both bounds
    b = _copy1d(a, astart, aend, 95, 110)
    assert np.array_equal(b, [0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0])


def test_parse_loci():
    loci = read_gtf('parse_loci.gtf')
    assert len(loci) == 3
    assert loci[0][0] == ('chrTest1', 10, 50)
    assert loci[1][0] == ('chrTest1', 50, 200)
    assert loci[2][0] == ('chrTest2', 100, 200)


def test_create_locus():
    loci = read_gtf('splice_sites.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 10, 525)
    t_dict = Transfrag.parse_gtf(gtf_lines)

    locus = Locus.create(t_dict.values())
    assert locus.chrom == 'chr1'
    assert locus.start == 10
    assert locus.end == 525
    a = locus.get_expr_data(49, 51, Strand.POS)
    assert np.array_equal(a, [1.0, 2.0])
    a = locus.get_expr_data(150, 151, Strand.POS)
    assert np.array_equal(a, [1.0])
    a = locus.get_expr_data(499, 501, Strand.POS)
    assert np.array_equal(a, [3.0, 1.0])
    with pytest.raises(TacoError):
        locus.get_expr_data(5, 15, Strand.POS)
    with pytest.raises(TacoError):
        locus.get_expr_data(520, 530, Strand.POS)
