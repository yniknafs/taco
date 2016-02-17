'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import networkx as nx

from taco.lib.base import Strand
from taco.lib.assemble import assemble_isoforms, Config
from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import create_optimal_path_graph, \
    create_path_graph
from taco.lib.path_finder import find_paths
import taco.lib.cpathfinder as cpathfinder

from taco.test.base import read_single_locus


def test_empty_graph_bug():
    t_dict, locus = read_single_locus('empty_graph_bug.gtf')
    transfrags = locus.get_transfrags(Strand.POS)
    sgraph = SpliceGraph.create(transfrags)
    isoforms = assemble_isoforms(sgraph, Config.defaults())
    assert len(isoforms) == 0


def test_path1():
    t_dict, locus = read_single_locus('path1.gtf')
    transfrags = locus.get_transfrags(Strand.POS)
    sgraph = SpliceGraph.create(transfrags)
    k = 2
    K = create_path_graph(sgraph, k)
    paths1 = find_paths(K, 'expr')
    paths2 = cpathfinder.find_paths(K, 'expr')
    assert len(paths1) == len(paths2)
    for p1, p2 in zip(paths1, paths2):
        p1, e1 = p1
        p2, e2 = p2
        assert p1 == p2
        assert abs(e1-e2) < 1e-8
    return


def test_path2():
    t_dict, locus = read_single_locus('noc2l_locus.gtf')
    for sgraph in locus.create_splice_graphs():
        K, k = create_optimal_path_graph(sgraph)
        paths1 = find_paths(K, 'expr')
        paths2 = cpathfinder.find_paths(K, 'expr')
        assert len(paths1) == len(paths2)
        for p1, p2 in zip(paths1, paths2):
            p1, e1 = p1
            p2, e2 = p2
            assert p1 == p2
            assert abs(e1-e2) < 1e-5


def test_topological_sort():
    t_dict, locus = read_single_locus('noc2l_locus.gtf')
    for sgraph in locus.create_splice_graphs():
        nodes1 = tuple(nx.topological_sort(sgraph.G))
        nodes2 = tuple(cpathfinder.topological_sort(sgraph.G))
        assert nodes1 == nodes2
