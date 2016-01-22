'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np

from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.transfrag import Transfrag
from taco.lib.splice_graph import aggregate_expression_data, \
    find_change_points, find_splice_sites, SpliceGraph, split_transfrag
from taco.lib.cChangePoint import find_change_points as c_find_change_points
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
        splice_sites.update(find_splice_sites(t))
    splice_sites = tuple(sorted(splice_sites))
    assert splice_sites == (100, 200, 250, 300, 400)
    # internal zero sites only
    expr_data = np.zeros((locus.end - locus.start), dtype=FLOAT_DTYPE)
    for t in transfrags:
        aggregate_expression_data(t, expr_data, locus.start)
    zero_sites = tuple(find_change_points(expr_data, locus.start))
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


def test_ref_starts_ends():
    t_dict, locus = read_single_locus('splice_graph.gtf')
    sg = SpliceGraph.create(t_dict.values())
    assert tuple(sorted(sg.ref_start_sites)) == (95,)
    assert tuple(sorted(sg.ref_stop_sites)) == (200,)
