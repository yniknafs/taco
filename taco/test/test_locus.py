'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import Strand, TacoError
from taco.lib.transfrag import Transfrag
from taco.lib.assemble import find_boundaries, split_transfrag, Locus
from taco.test.base import read_gtf


def test_parse_loci():
    loci = read_gtf('parse_loci.gtf')
    assert len(loci) == 3
    assert loci[0][0] == ('chrTest1', 10, 50)
    assert loci[1][0] == ('chrTest1', 50, 200)
    assert loci[2][0] == ('chrTest2', 100, 200)


def test_find_splice_sites():
    loci = read_gtf('splice_sites.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 10, 525)
    t_dict = Transfrag.parse_gtf(gtf_lines)
    splice_sites = tuple(find_boundaries(t_dict.values()))
    assert splice_sites == (10, 100, 200, 250, 300, 400, 525)


def test_find_transfrag_nodes():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    boundaries = find_boundaries(t_dict.values())
    # check nodes
    t = t_dict['A']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (200, 250), (250, 300), (400, 525))
    t = t_dict['B']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (250, 300), (400, 525))
    t = t_dict['C']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((100, 200), (200, 250), (250, 300), (400, 525))
    t = t_dict['D']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((300, 400), (400, 525))


def test_create_locus():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    locus = Locus.create(t_dict.values())
    assert tuple(locus.boundaries) == (10, 100, 200, 250, 300, 400, 525)
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
