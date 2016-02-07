'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import networkx as nx
import matplotlib.pyplot as plt

from taco.lib.base import Exon, Strand
from taco.lib.transfrag import Transfrag

from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import choose_k_by_frag_length, create_path_graph, \
    get_kmers, add_path, SOURCE, SINK, smooth_graph, get_unreachable_kmers, \
    get_path, create_optimal_path_graph

from taco.test.base import read_single_locus

def test_cluster():
    pass


def cluster_isoforms(paths):

    node_cluster_map = {}
    clusters = {}
    # initialize each node to be in a "cluster" by itself
    for n in G.nodes_iter():
        node_chain_map[n] = n
        chains[n] = set((n,))
    for u,v in G.edges_iter():
        if not can_collapse_func(G,u,v):
            continue
        # get chains containing these nodes
        u_new = node_chain_map[u]
        u_chain = chains[u_new]
        del chains[u_new]
        v_new = node_chain_map[v]
        v_chain = chains[v_new]
        del chains[v_new]
        # merge chains
        merged_chain = u_chain.union(v_chain)
        merged_node = Exon(imin2(u_new.start, v_new.start),
                           imax2(u_new.end, v_new.end))
        # point all nodes in chain to new parent
        for n in merged_chain:
            node_chain_map[n] = merged_node
        chains[merged_node] = merged_chain
    # sort chain nodes by genome position and store as list
    for parent in chains:
        chains[parent] = sorted(chains[parent], key=operator.attrgetter('start'))
    return node_chain_map, chains
