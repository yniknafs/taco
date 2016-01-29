'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import networkx as nx
import matplotlib.pyplot as plt

from taco.lib.base import Exon, Strand
from taco.lib.transfrag import Transfrag

from taco.lib.splice_graph import find_node_boundaries
from taco.lib.path_graph import choose_k
from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import create_path_graph, get_kmers, add_path, \
    SOURCE, SINK

from taco.test.base import read_single_locus


def testchoose_k():
    A = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 100), Exon(200, 300), Exon(400, 500)],
                  _id='A')
    node_bounds = find_node_boundaries([A])
    # case does not reach min frag length
    k = choose_k([A], node_bounds, min_path_length=400)
    assert k == 3
    # case reaches min frag length
    k = choose_k([A], node_bounds, min_path_length=300)
    assert k == 3
    k = choose_k([A], node_bounds, min_path_length=250)
    assert k == 3
    # case reaches frag length with nodes to spare
    k = choose_k([A], node_bounds, min_path_length=200)
    assert k == 2
    # case reaches frag length with nodes to spare
    k = choose_k([A], node_bounds, min_path_length=100)
    assert k == 1


def testchoose_k_2():
    A = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10),
                         Exon(20, 30),
                         Exon(40, 50),
                         Exon(60, 70),
                         Exon(80, 90),
                         Exon(100, 110),
                         Exon(110, 120),
                         Exon(130, 140),
                         Exon(150, 160),
                         Exon(170, 180),
                         Exon(190, 200)],
                  _id='A')
    B = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10),
                         Exon(500, 600),
                         Exon(1000, 1100)],
                  _id='B')
    node_bounds = find_node_boundaries([A, B])
    k = choose_k([A, B], node_bounds, min_path_length=200)
    assert k == 11
    k = choose_k([A, B], node_bounds, min_path_length=100)
    assert k == 10
    k = choose_k([A, B], node_bounds, min_path_length=1)
    assert k == 1

    B = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10),
                         Exon(65, 75),
                         Exon(500, 600),
                         Exon(1000, 1100)],
                  _id='B')
    node_bounds = find_node_boundaries([A, B])
    k = choose_k([A, B], node_bounds, min_path_length=200)
    assert k == 12
    k = choose_k([A, B], node_bounds, min_path_length=100)
    assert k == 11
    k = choose_k([A, B], node_bounds, min_path_length=50)
    assert k == 6
    k = choose_k([A, B], node_bounds, min_path_length=1)
    assert k == 1


def test_path_graph1():
    # read transcripts
    t_dict, locus = read_single_locus('path1.gtf')
    SG = SpliceGraph.create(t_dict.values())
    # paths
    ABCDE = (SOURCE, Exon(0, 100), Exon(200, 300), Exon(400, 500),
             Exon(600, 700), Exon(800, 900), SINK)
    ACE = (SOURCE, Exon(0, 100), Exon(400, 500), Exon(800, 900), SINK)
    ABCE = (SOURCE, Exon(0, 100), Exon(200, 300), Exon(400, 500),
            Exon(800, 900), SINK)
    ACDE = (SOURCE, Exon(0, 100), Exon(400, 500), Exon(600, 700),
            Exon(800, 900), SINK)
    paths = [ABCDE, ACE, ABCE, ACDE]
    # create path graph k = 2
    k = 2
    G1, lost_paths = create_path_graph(SG, k)
    G2 = nx.DiGraph()
    for path in paths:
        kmers = list(get_kmers(path, k))
        add_path(G2, kmers, 1.0)
    assert nx.is_isomorphic(G1, G2)
