'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.transfrag import Transfrag
from taco.lib.splice_graph import aggregate_expression_data, find_zero_sites, \
    find_splice_sites, SpliceGraph, split_transfrag
from taco.lib.base import Strand
from taco.test.base import read_gtf, read_single_locus


def test_find_zero_sites():
    a = np.ones(10)
    assert len(find_zero_sites(a)) == 0
    a = np.zeros(10)
    with pytest.raises(AssertionError):
        find_zero_sites(a)
    a = np.array([0, 1, 1, 0])
    with pytest.raises(AssertionError):
        find_zero_sites(a)
    a = np.array([1, 0, 0, 0, 1])
    assert tuple(find_zero_sites(a)) == (1, 4)
    a = np.array([1, 0, 1, 0, 1])
    assert tuple(find_zero_sites(a)) == (1, 2, 3, 4)


def test_find_boundaries():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    chrom, start, end = interval
    transfrags = Transfrag.parse_gtf(gtf_lines).values()
    # splice sites only
    splice_sites = tuple(find_splice_sites(transfrags))
    assert splice_sites == (100, 200, 250, 300, 400)
    # internal zero sites only
    a = aggregate_expression_data(transfrags, start, end)
    zero_sites = tuple(find_zero_sites(a, start))
    assert zero_sites == (100, 150, 300, 375)
    # combined boundaries
    sg = SpliceGraph.create(transfrags)
    assert tuple(sg.boundaries) == (10, 100, 150, 200, 250, 300, 375, 400, 525)


def test_split_transfrag():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    sg = SpliceGraph.create(t_dict.values())
    boundaries = sg.boundaries
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
