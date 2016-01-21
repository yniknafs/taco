'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np

from taco.lib.base import Strand
from taco.test.base import read_single_locus


def test_impute_strand():
    t_dict, locus = read_single_locus('impute_strand.gtf')
    assert len(locus.strand_transfrags[Strand.POS]) == 1
    assert len(locus.strand_transfrags[Strand.NEG]) == 1
    assert len(locus.strand_transfrags[Strand.NA]) == 3
    locus.impute_unknown_strands()
    assert len(locus.strand_transfrags[Strand.POS]) == 2
    assert len(locus.strand_transfrags[Strand.NEG]) == 2
    assert len(locus.strand_transfrags[Strand.NA]) == 1
    a = locus.get_expr_data(9, 11, Strand.POS)
    assert np.array_equal(a, [2.0, 1.0])
    a = locus.get_expr_data(14, 16, Strand.POS)
    assert np.array_equal(a, [1.0, 0.0])
    a = locus.get_expr_data(14, 16, Strand.NEG)
    assert np.array_equal(a, [0.0, 1.0])
    a = locus.get_expr_data(14, 16, Strand.NA)
    assert np.array_equal(a, [1.0, 1.0])
    a = locus.get_expr_data(19, 21, Strand.NEG)
    assert np.array_equal(a, [1.0, 2.0])


def test_impute_strand_guided():
    t_dict, locus = read_single_locus('impute_strand_guided.gtf',
                                      guided_strand=True)
    assert len(locus.strand_transfrags[Strand.POS]) == 2
    assert len(locus.strand_transfrags[Strand.NEG]) == 1
    assert len(locus.strand_transfrags[Strand.NA]) == 3
    locus.impute_unknown_strands()
    assert t_dict['C'].strand == Strand.POS
    assert len(locus.strand_transfrags[Strand.POS]) == 4
    assert len(locus.strand_transfrags[Strand.NEG]) == 1
    assert len(locus.strand_transfrags[Strand.NA]) == 1
    a = locus.get_expr_data(14, 16, Strand.POS)
    assert np.array_equal(a, [2.0, 1.0])
    a = locus.get_expr_data(14, 16, Strand.NEG)
    assert np.array_equal(a, [0.0, 0.0])
    a = locus.get_expr_data(14, 16, Strand.NA)
    assert np.array_equal(a, [0.0, 1.0])
    a = locus.get_expr_data(19, 21, Strand.NEG)
    assert np.array_equal(a, [0.0, 1.0])
    a = locus.get_expr_data(19, 21, Strand.NA)
    assert np.array_equal(a, [1.0, 1.0])
