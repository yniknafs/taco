'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from taco.lib.base import Strand
from taco.lib.assemble import assemble_isoforms, Config
from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_finder import find_paths
from taco.lib.path_graph import PathGraphFactory, PathGraph
from taco.lib.cpathfinder import find_paths

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
    pgf = PathGraphFactory(sgraph)
    pgraph = pgf.create(k)
    paths = find_paths(pgraph)
    return


def test_path2():
    t_dict, locus = read_single_locus('noc2l_locus.gtf')
    for sgraph in locus.create_splice_graphs():
        pgraphfactory = PathGraphFactory(sgraph)
        pgraph, k = pgraphfactory.create_optimal()
        paths = find_paths(pgraph)
    return


# def test_path_ties():
#     G = PathGraph()
#     G.add_path((G.SOURCE, 20, 30, 40, 50, 60, 70, 80, G.SINK), 10.0)
#     G.add_path((50, 70), 50.0)
#     paths = find_paths2(G)
#     assert len(paths) == 1
#     p, e = paths[0]
#     p = tuple(G.nodes[i] for i in p)
#     assert p == (-1, 20, 30, 40, 50, 70, 80, -2)
#     assert e == 10.0
