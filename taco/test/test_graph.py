'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from taco.lib.graph import Graph
from taco.lib.path_graph2 import PathGraphFactory

from taco.test.base import read_single_locus


def test_add_paths():
    G = Graph()
    G.add_path((Graph.SOURCE, 10, 20, 30, 40, Graph.SINK))
    G.add_path((Graph.SOURCE, 10, 30, 40, Graph.SINK))
    G.add_path((Graph.SOURCE, 10, Graph.SINK))
    assert len(G.preds[Graph.SOURCE_ID]) == 0
    assert len(G.succs[Graph.SINK_ID]) == 0
    assert sorted(G.nodes[x] for x in G.succs[Graph.SOURCE_ID]) == [10]
    assert sorted(G.nodes[x] for x in G.preds[Graph.SINK_ID]) == [10, 40]


def test_remove_node():
    G = Graph()
    G.add_path((Graph.SOURCE, 10, 20, 30, 40, Graph.SINK))
    G.add_path((Graph.SOURCE, 10, 30, 40, Graph.SINK))
    G.add_path((Graph.SOURCE, 10, Graph.SINK))
    assert len(G) == 4
    assert len(G.preds[G.node_id_map[20]]) == 1
    assert len(G.preds[G.node_id_map[30]]) == 2
    node_id = G.node_id_map[10]
    G.remove_node_id(G.node_id_map[10])
    assert 10 not in G.node_id_map
    assert G.nodes[node_id] == Graph.EMPTY
    assert len(G) == 3
    assert len(G.preds[G.node_id_map[20]]) == 0
    assert len(G.preds[G.node_id_map[30]]) == 1


def test_unreachable():
    G = Graph()
    G.add_path((Graph.SOURCE,) + tuple(range(10)) + (Graph.SINK,))
    assert len(G.get_unreachable_nodes()) == 0
    # unreachable from sink
    G.add_path((2, 20, 21))
    unreachable = set(G.nodes[x] for x in G.get_unreachable_nodes())
    assert unreachable == set([20, 21])
    # unreachable from source
    G.add_path((33, 32, 31, 9))
    unreachable = set(G.nodes[x] for x in G.get_unreachable_nodes())
    assert unreachable == set([20, 21, 33, 32, 31])


def test_is_valid():
    G = Graph()
    assert not G.is_valid()
    G.add_path((Graph.SOURCE,) + tuple(range(10)) + (Graph.SINK,))
    assert G.is_valid()
    G.remove_node_id(G.node_id_map[5])
    assert not G.is_valid()
    G.remove_unreachable_nodes()
    assert not G.is_valid()


def test_is_topological_sort():
    G = Graph()
    G.add_path((G.SOURCE, 10, 20, 30, 40, G.SINK))
    G.add_path((G.SOURCE, 10, 30, 40, G.SINK))
    G.add_path((G.SOURCE, 10, G.SINK))
    G.add_path((G.SOURCE, 20, G.SINK))

    toposort = (G.SOURCE, 10, 30, 20, 40, G.SINK)
    order = (G.node_id_map[x] for x in toposort)
    assert not G.is_topological_sort(order)
    toposort = (G.SOURCE, 20, 10, 30, 40, G.SINK)
    order = (G.node_id_map[x] for x in toposort)
    assert not G.is_topological_sort(order)
    toposort = (G.SOURCE, 10, 20, 30, 40, G.SINK)
    order = (G.node_id_map[x] for x in toposort)
    assert G.is_topological_sort(order)

    G = Graph()
    G.add_path((G.SOURCE, 10, 30, 40, G.SINK))
    G.add_path((G.SOURCE, 10, 20, 40, G.SINK))

    toposort = (G.SOURCE, 10, 20, 30, 40, G.SINK)
    order = (G.node_id_map[x] for x in toposort)
    assert G.is_topological_sort(order)
    toposort = (G.SOURCE, 10, 30, 20, 40, G.SINK)
    order = (G.node_id_map[x] for x in toposort)
    assert G.is_topological_sort(order)


def test_topological_sort():
    G = Graph()
    G.add_path((G.SOURCE, 10, 20, 30, 40, G.SINK))
    G.add_path((G.SOURCE, 10, 30, 40, G.SINK))
    G.add_path((G.SOURCE, 10, G.SINK))
    G.add_path((G.SOURCE, 20, G.SINK))
    assert G.is_topological_sort(G.topological_sort())
    t_dict, locus = read_single_locus('noc2l_locus.gtf')
    for sgraph in locus.create_splice_graphs():
        pgf = PathGraphFactory(sgraph)
        G = pgf.create(k=1)
        assert G.is_topological_sort(G.topological_sort())
        assert G.is_topological_sort(G.topological_sort_dfs())
