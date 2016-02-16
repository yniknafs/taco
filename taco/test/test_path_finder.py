'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

from taco.lib.base import Strand
from taco.lib.assemble import assemble_isoforms, Config
from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import create_optimal_path_graph
from taco.lib.path_finder import find_paths

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
    K, k = create_optimal_path_graph(sgraph)
    paths1 = list(find_paths(K, K.graph['source'], K.graph['sink']))
