'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import TacoError, Strand, Exon
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.transfrag import Transfrag
from taco.lib.cchangepoint import find_threshold_points as \
    find_threshold_points_cy
from taco.lib.changepoint import find_threshold_points
from taco.lib.locus import Locus
from taco.lib.splice_graph import split_transfrag, SpliceGraph, SGNode

from taco.test.base import read_gtf, read_single_locus


def test_find_threshold_points():
    for func in (find_threshold_points, find_threshold_points_cy):
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


def test_find_node_boundaries():
    t_dict, locus = read_single_locus('splice_sites.gtf')
    transfrags = t_dict.values()
    splice_sites = set()
    for t in transfrags:
        splice_sites.update(t.itersplices())
    splice_sites = tuple(sorted(splice_sites))
    assert splice_sites == (100, 200, 250, 300, 400)
    # aggregate expression
    sg = SpliceGraph.create(transfrags)
    # zero change points
    zero_sites = tuple(find_threshold_points(sg.expr_data, sg.start))
    assert zero_sites == (100, 150, 300, 375)
    # combined boundaries
    boundaries = tuple(sg._find_node_boundaries())
    assert boundaries == (10, 100, 150, 200, 250, 300, 375, 400, 525)


def test_ref_starts_ends():
    t_dict, locus = read_single_locus('change_point1.gtf')
    sg = SpliceGraph.create(t_dict.values())
    assert tuple(sorted(sg.ref_start_sites)) == (95,)
    assert tuple(sorted(sg.ref_stop_sites)) == (200,)


def test_mark_start_stop_sites1():
    t_dict, locus = read_single_locus('change_point1.gtf')
    sgraph = SpliceGraph.create(t_dict.values())
    G = sgraph.G
    assert len(G) == 1
    n_id = sgraph.get_node_id(Exon(50, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    # add a start site change point
    sgraph.start_sites.add(125)
    sgraph.recreate()
    G = sgraph.G
    assert len(G) == 2
    n_id = sgraph.get_node_id(Exon(50, 125))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and not nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(125, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    # add a stop site change point
    sgraph.stop_sites.add(80)
    sgraph.recreate()
    G = sgraph.G
    assert len(G) == 3
    n_id = sgraph.get_node_id(Exon(50, 80))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(80, 125))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert not nd[SGNode.IS_START] and not nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(125, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]

    # flip strand
    for t_id, t in t_dict.iteritems():
        t.strand = Strand.NEG
    sgraph = SpliceGraph.create(t_dict.values())
    G = sgraph.G
    assert len(G) == 1

    n_id = sgraph.get_node_id(Exon(50, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    # add a start site change point
    sgraph.start_sites.add(125)
    sgraph.recreate()
    G = sgraph.G
    assert len(G) == 2
    n_id = sgraph.get_node_id(Exon(50, 125))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(125, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and not nd[SGNode.IS_STOP]

    # add a stop site change point
    sgraph.stop_sites.add(80)
    sgraph.recreate()
    G = sgraph.G
    assert len(G) == 3
    n_id = sgraph.get_node_id(Exon(50, 80))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert not nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(80, 125))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and nd[SGNode.IS_STOP]
    n_id = sgraph.get_node_id(Exon(125, 200))
    assert n_id in G
    nd = sgraph.G.node[n_id]
    assert nd[SGNode.IS_START] and not nd[SGNode.IS_STOP]


def test_mark_start_stop_sites2():
    # pos strand not guided
    t_dict, locus = read_single_locus('multi_strand1.gtf')
    sgraph = SpliceGraph.create(locus.get_transfrags(Strand.POS))
    G = sgraph.G
    assert G.node[sgraph.get_node_id(Exon(100, 200))][SGNode.IS_START]
    assert G.node[sgraph.get_node_id(Exon(400, 650))][SGNode.IS_STOP]

    # neg strand not guided
    sgraph = SpliceGraph.create(locus.get_transfrags(Strand.NEG))
    G = sgraph.G
    assert G.node[sgraph.get_node_id(Exon(950, 980))][SGNode.IS_START]
    assert G.node[sgraph.get_node_id(Exon(400, 500))][SGNode.IS_STOP]

    # pos strand guided
    sgraph = SpliceGraph.create(locus.get_transfrags(Strand.POS),
                                guided_ends=True,
                                guided_assembly=True)
    G = sgraph.G
    assert G.node[sgraph.get_node_id(Exon(100, 150))][SGNode.IS_START]
    assert G.node[sgraph.get_node_id(Exon(150, 200))][SGNode.IS_START]
    assert G.node[sgraph.get_node_id(Exon(500, 600))][SGNode.IS_STOP]
    assert G.node[sgraph.get_node_id(Exon(600, 650))][SGNode.IS_STOP]
    assert G.node[sgraph.get_node_id(Exon(150, 200))][SGNode.IS_REF]
    assert G.node[sgraph.get_node_id(Exon(300, 400))][SGNode.IS_REF]
    assert G.node[sgraph.get_node_id(Exon(500, 600))][SGNode.IS_REF]
    assert not G.node[sgraph.get_node_id(Exon(100, 150))][SGNode.IS_REF]
    assert not G.node[sgraph.get_node_id(Exon(600, 650))][SGNode.IS_REF]

    # neg strand guided
    sgraph = SpliceGraph.create(locus.get_transfrags(Strand.NEG),
                                guided_ends=True,
                                guided_assembly=True)
    G = sgraph.G
    assert G.node[sgraph.get_node_id(Exon(350, 400))][SGNode.IS_STOP]
    assert G.node[sgraph.get_node_id(Exon(980, 1000))][SGNode.IS_START]
    assert not G.node[sgraph.get_node_id(Exon(950, 980))][SGNode.IS_START]
    for n, nd in G.nodes_iter(data=True):
        assert nd[SGNode.IS_REF]
    return


def test_split_transfrag():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    sg = SpliceGraph.create(t_dict.values())
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


def test_multi_strand1():
    # read gtf and test basic values
    loci = read_gtf('multi_strand1.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 100, 1000)
    t_dict = Transfrag.parse_gtf(gtf_lines)
    assert len(t_dict) == 5
    locus = Locus.create(t_dict.values())
    assert locus.chrom == 'chr1'
    assert locus.start == 100
    assert locus.end == 1000
    # raise exception when creating with multiple strands
    with pytest.raises(TacoError):
        SpliceGraph.create(t_dict.values())
    transfrags_pos = locus.get_transfrags(Strand.POS)
    transfrags_neg = locus.get_transfrags(Strand.NEG)
    sgpos = SpliceGraph.create(transfrags_pos)
    sgneg = SpliceGraph.create(transfrags_neg)

    # test
    assert sgpos.chrom == 'chr1'
    assert sgpos.start == 100
    assert sgpos.end == 650
    assert sgpos.strand == Strand.POS
    assert sgpos.ref_start_sites == [150]
    assert sgpos.ref_stop_sites == [600]
    with pytest.raises(TacoError):
        sgpos.get_expr_data(90, 110)
    with pytest.raises(TacoError):
        sgpos.get_expr_data(650, 655)
    assert np.array_equal(sgpos.get_expr_data(100, 105), np.ones(5))

    assert sgneg.chrom == 'chr1'
    assert sgneg.start == 350
    assert sgneg.end == 1000
    assert sgneg.strand == Strand.NEG
    assert sgneg.ref_start_sites == [1000]
    assert sgneg.ref_stop_sites == [350]
    with pytest.raises(TacoError):
        sgneg.get_expr_data(340, 350)
    with pytest.raises(TacoError):
        sgneg.get_expr_data(1000, 1010)
    assert np.array_equal(sgneg.get_expr_data(400, 405), np.ones(5))
    assert np.array_equal(sgneg.get_expr_data(945, 950), np.zeros(5))
    assert np.array_equal(sgneg.get_expr_data(950, 955), np.ones(5))
    assert np.array_equal(sgneg.get_expr_data(980, 985), np.zeros(5))

    # test locus boundaries
    bpos = tuple(sgpos._find_node_boundaries())
    assert bpos == tuple((100, 200, 300, 400, 650))
    bneg = tuple(sgneg._find_node_boundaries())
    assert bneg == tuple((350, 400, 500, 950, 980, 1000))

    # added guided ends/assembly to use boundaries from reference
    lpos = SpliceGraph.create(transfrags_pos,
                              guided_ends=True,
                              guided_assembly=True)
    bpos = tuple(lpos._find_node_boundaries())
    assert bpos == tuple((100, 150, 200, 300, 400, 500, 600, 650))

    lneg = SpliceGraph.create(transfrags_neg,
                              guided_ends=True,
                              guided_assembly=True)
    bneg = tuple(lneg._find_node_boundaries())
    assert bneg == tuple((350, 400, 500, 750, 900, 950, 980, 1000))


def test_multi_strand2():
    t_dict, locus = read_single_locus('multi_strand2.gtf')
    transfrags_pos = locus.get_transfrags(Strand.POS)
    sgpos = SpliceGraph.create(transfrags_pos)
    sgdict = {}
    for sg in sgpos.split():
        k = ('%s:%d-%d[%s]' % (sg.chrom, sg.start, sg.end,
             Strand.to_gtf(sg.strand)))
        sgdict[k] = sg
    assert 'chr1:100-300[+]' in sgdict
    assert 'chr1:400-600[+]' in sgdict
