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


def test_reachability():
    source = -1
    sink = -2
    G = nx.DiGraph()
    nx.path_graph(10, create_using=G)
    G.add_edge(source, 0)
    G.add_edge(9, sink)
    assert len(get_unreachable_kmers(G, source, sink)) == 0
    # unreachable from sink
    G.add_edges_from([(2, 20), (20, 21)])
    unreachable = get_unreachable_kmers(G, source, sink)
    assert unreachable == set([20, 21])
    # unreachable from source
    G.add_edges_from([(33, 32), (32, 31), (31, 9)])
    unreachable = get_unreachable_kmers(G, source, sink)
    assert unreachable == set([20, 21, 33, 32, 31])


def test_unreachable_kmers():
    t_dict, locus = read_single_locus('path_graph_k2.gtf')
    sgraph = SpliceGraph.create(t_dict.values())
    K = create_path_graph(sgraph, k=2)
    assert not K.graph['valid']
    assert len(K) == 0
    assert len(K.graph['unreachable_kmers']) == 8

    K = create_path_graph(sgraph, k=1)
    assert len(K.graph['lost_transfrags']) == 0
    assert len(K.graph['unreachable_kmers']) == 0
    assert K.graph['valid']
    assert len(K) == 8

    K, k = create_optimal_path_graph(sgraph, frag_length=0, kmax=0,
                                     loss_threshold=1.0)
    assert k == 1
    assert len(K) == 8

    # print 'K', K.nodes()
    # print 'lost', K.graph['lost_transfrags']
    # print 'unreachable', unreachable
    # print 'valid', valid
    #
    # print 'start stop', sgraph.get_start_stop_nodes()
    # for t in sgraph.itertransfrags():
    #     print t._id, t.is_ref, t.expr, t.exons, get_path(sgraph, t)
    #
    # for n in K:
    #     if n < 0: continue
    #     print 'node', n, 'kmer', K.graph['id_kmer_map'][n], K.node[n]['expr']
    #
    # nx.draw_networkx(K, with_labels=True)
    # plt.show()


def test_choose_k_1():
    A = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 100), Exon(200, 300), Exon(400, 500)],
                  _id='A')
    sgraph = SpliceGraph.create([A])
    # case does not reach min frag length
    k = choose_k_by_frag_length(sgraph, 400, 1)
    assert k == 3
    # case reaches min frag length
    k = choose_k_by_frag_length(sgraph, 300, 1)
    assert k == 3
    k = choose_k_by_frag_length(sgraph, 250, 1)
    assert k == 3
    # case reaches frag length with nodes to spare
    k = choose_k_by_frag_length(sgraph, 200, 1)
    assert k == 2
    # case reaches frag length with nodes to spare
    k = choose_k_by_frag_length(sgraph, 100, 1)
    assert k == 1
    # enforce kmin
    k = choose_k_by_frag_length(sgraph, 100, kmin=2)
    assert k == 2
    k = choose_k_by_frag_length(sgraph, 100, kmin=3)
    assert k == 3
    k = choose_k_by_frag_length(sgraph, 100, kmin=5)
    assert k == 3


def test_choose_k_2():
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
    sgraph = SpliceGraph.create([A, B])
    k = choose_k_by_frag_length(sgraph, frag_length=200, kmin=1)
    assert k == 11
    k = choose_k_by_frag_length(sgraph, frag_length=100, kmin=1)
    assert k == 10
    k = choose_k_by_frag_length(sgraph, frag_length=1, kmin=1)
    assert k == 1
    k = choose_k_by_frag_length(sgraph, frag_length=1, kmin=3)
    assert k == 3

    B = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10),
                         Exon(65, 75),
                         Exon(500, 600),
                         Exon(1000, 1100)],
                  _id='B')
    sgraph = SpliceGraph.create([A, B])
    k = choose_k_by_frag_length(sgraph, frag_length=200)
    assert k == 12
    k = choose_k_by_frag_length(sgraph, frag_length=100)
    assert k == 11
    k = choose_k_by_frag_length(sgraph, frag_length=50)
    assert k == 6
    k = choose_k_by_frag_length(sgraph, frag_length=1, kmin=1)
    assert k == 1
    k = choose_k_by_frag_length(sgraph, frag_length=1, kmin=2)
    assert k == 2


def test_choose_k_3():
    A = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10000),
                         Exon(20000, 30000),
                         Exon(40000, 50000),
                         Exon(60000, 70000),
                         Exon(80000, 90000),
                         Exon(100000, 110000),
                         Exon(120000, 130000),
                         Exon(140000, 150000),
                         Exon(160000, 170000),
                         Exon(180000, 190000)],
                  _id='A')
    sgraph = SpliceGraph.create([A])
    k = choose_k_by_frag_length(sgraph, frag_length=400, kmin=1)
    assert k == 1
    k = choose_k_by_frag_length(sgraph, frag_length=400, kmin=2)
    assert k == 2
    k = choose_k_by_frag_length(sgraph, frag_length=20000, kmin=2)
    assert k == 2
    k = choose_k_by_frag_length(sgraph, frag_length=100000, kmin=2)
    assert k == 10


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
    G1 = create_path_graph(SG, k)
    G2 = nx.DiGraph()
    for path in paths:
        kmers = list(get_kmers(path, k))
        add_path(G2, kmers, 1.0)
    assert nx.is_isomorphic(G1, G2)


def test_path_graph2():
    t_dict, locus = read_single_locus('change_point2.gtf')
    sgraph = SpliceGraph.create(t_dict.values())

    # trivial case without additional stops or starts
    k = choose_k_by_frag_length(sgraph, frag_length=400, kmin=2)
    assert k == 1
    K = create_path_graph(sgraph, k)
    kmer_id_map = dict((v, k) for k, v in K.graph['id_kmer_map'].iteritems())
    assert len(K.graph['lost_transfrags']) == 0
    n = kmer_id_map[(Exon(0, 100),)]
    assert K.node[n]['expr'] == 12.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0

    # add a stop site
    sgraph.stop_sites.add(50)
    sgraph.recreate()
    K = create_path_graph(sgraph, k=2)
    assert len(K.graph['lost_transfrags']) == 0
    kmer_id_map = dict((v, k) for k, v in K.graph['id_kmer_map'].iteritems())
    n1 = kmer_id_map[(Exon(start=0, end=50), Exon(start=50, end=100))]
    n2 = kmer_id_map[(Exon(start=0, end=50),)]
    assert K.node[n1]['expr'] == 2.0
    assert K.node[n2]['expr'] == 10.0
    assert K.node[SOURCE]['expr'] == 11.0
    assert K.node[SINK]['expr'] == 11.0
    # smooth kmer graph
    smooth_graph(K)
    assert K.node[n1]['expr'] == 2.0
    assert K.node[n2]['expr'] == 10.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0

    # add both a start and a stop site
    sgraph.start_sites.add(50)
    sgraph.stop_sites.add(50)
    sgraph.recreate()
    K = create_path_graph(sgraph, k=2)
    smooth_graph(K)
    kmer_id_map = dict((v, k) for k, v in K.graph['id_kmer_map'].iteritems())
    n1 = kmer_id_map[(Exon(start=0, end=50), Exon(start=50, end=100))]
    n2 = kmer_id_map[(Exon(start=0, end=50),)]
    n3 = kmer_id_map[(Exon(start=50, end=100),)]
    assert K.node[n1]['expr'] == 1.0
    assert K.node[n2]['expr'] == 10.0
    assert K.node[n3]['expr'] == 1.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0
