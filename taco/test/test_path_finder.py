'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import networkx as nx

from taco.lib.path_graph import init_node_attrs, SOURCE, SINK
from taco.lib.path_finder import find_path, init_tmp_attributes


def test_path_finder1():
    G = nx.DiGraph()
    # nodes
    for n in [SOURCE, 1, 2, SINK]:
        G.add_node(n, attr_dict=init_node_attrs())
    # edges
    G.add_edges_from([(SOURCE, 1), (1, 2), (2, SINK)])
    # find path
    init_tmp_attributes(G)
    print find_path(G, SOURCE, SINK)
