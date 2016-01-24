'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import TacoError, Strand
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.transfrag import Transfrag
from taco.lib.cChangePoint import find_change_points as c_find_change_points
from taco.lib.changepoint import find_change_points
from taco.lib.locus import Locus, StrandedLocus, split_transfrag

from taco.test.base import read_gtf, read_single_locus


def test_find_change_points():
    for func in (find_change_points, c_find_change_points):
        a = np.ones(10, dtype=FLOAT_DTYPE)
        assert len(func(a)) == 0
        a = np.zeros(10, dtype=FLOAT_DTYPE)
        assert len(func(a)) == 0
        a = np.array([0, 1, 1, 0], dtype=FLOAT_DTYPE)
        assert tuple(func(a)) == (1, 3)
        a = np.array([1, 0, 0, 0, 1], dtype=FLOAT_DTYPE)
        assert tuple(func(a)) == (1, 4)
        a = np.array([1, 0, 1, 0, 1], dtype=FLOAT_DTYPE)
        assert tuple(func(a)) == (1, 2, 3, 4)


def test_find_boundaries():
    t_dict, locus = read_single_locus('splice_sites.gtf')
    transfrags = t_dict.values()
    splice_sites = set()
    for t in transfrags:
        splice_sites.update(t.itersplices())
    splice_sites = tuple(sorted(splice_sites))
    assert splice_sites == (100, 200, 250, 300, 400)
    # aggregate expression
    slocus = StrandedLocus.create(transfrags)
    # zero change points
    zero_sites = tuple(find_change_points(slocus.expr_data, slocus.start))
    assert zero_sites == (100, 150, 300, 375)
    # combined boundaries
    boundaries = tuple(slocus._find_node_boundaries())
    assert boundaries == (10, 100, 150, 200, 250, 300, 375, 400, 525)


def test_split_transfrag():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    sg = StrandedLocus.create(t_dict.values())
    boundaries = tuple(sg._find_node_boundaries())
    # check nodes
    t = t_dict['A']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (200, 250), (250, 300), (400, 525))
    t = t_dict['B']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (250, 300), (400, 525))
    t = t_dict['C']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((150, 200), (200, 250), (250, 300), (400, 525))
    t = t_dict['D']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((375, 400), (400, 525))


def test_multi_strand():
    # read gtf and test basic values
    loci = read_gtf('locus_multi_strand.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 100, 1000)
    t_dict = Transfrag.parse_gtf(gtf_lines)
    assert len(t_dict) == 5
    locus = Locus.create(t_dict.values())
    assert locus.chrom == 'chr1'
    assert locus.start == 100
    assert locus.end == 1000
    # raise exception when creating StrandedLocus with multiple strands
    with pytest.raises(TacoError):
        StrandedLocus.create(t_dict.values())
    transfrags_pos = locus.get_transfrags(Strand.POS)
    transfrags_neg = locus.get_transfrags(Strand.NEG)
    # create StrandedLocus
    slocus_pos = StrandedLocus.create(transfrags_pos)
    slocus_neg = StrandedLocus.create(transfrags_neg)

    # test locus boundaries and expression levels
    assert slocus_pos.chrom == 'chr1'
    assert slocus_pos.start == 100
    assert slocus_pos.end == 650
    assert slocus_pos.strand == Strand.POS
    assert slocus_pos.ref_start_sites == [150]
    assert slocus_pos.ref_stop_sites == [600]
    with pytest.raises(AssertionError):
        slocus_pos.get_expr_data(90, 110)
    with pytest.raises(AssertionError):
        slocus_pos.get_expr_data(650, 655)
    assert np.array_equal(slocus_pos.get_expr_data(100, 105), np.ones(5))

    assert slocus_neg.chrom == 'chr1'
    assert slocus_neg.start == 350
    assert slocus_neg.end == 1000
    assert slocus_neg.strand == Strand.NEG
    assert slocus_neg.ref_start_sites == [1000]
    assert slocus_neg.ref_stop_sites == [350]
    with pytest.raises(AssertionError):
        slocus_neg.get_expr_data(340, 350)
    with pytest.raises(AssertionError):
        slocus_neg.get_expr_data(1000, 1010)
    assert np.array_equal(slocus_neg.get_expr_data(400, 405), np.ones(5))
    assert np.array_equal(slocus_neg.get_expr_data(945, 950), np.zeros(5))
    assert np.array_equal(slocus_neg.get_expr_data(950, 955), np.ones(5))
    assert np.array_equal(slocus_neg.get_expr_data(980, 985), np.zeros(5))

    # test locus boundaries
    bpos = tuple(slocus_pos._find_node_boundaries())
    assert bpos == tuple((100, 200, 300, 400, 650))
    bneg = tuple(slocus_neg._find_node_boundaries())
    assert bneg == tuple((350, 400, 500, 950, 980, 1000))

    # added guided ends/assembly to use boundaries from reference
    lpos = StrandedLocus.create(transfrags_pos,
                                guided_ends=True,
                                guided_assembly=True)
    bpos = tuple(lpos._find_node_boundaries())
    assert bpos == tuple((100, 150, 200, 300, 400, 500, 600, 650))

    lneg = StrandedLocus.create(transfrags_neg,
                                guided_ends=True,
                                guided_assembly=True)
    bneg = tuple(lneg._find_node_boundaries())
    assert bneg == tuple((350, 400, 500, 750, 900, 950, 980, 1000))
