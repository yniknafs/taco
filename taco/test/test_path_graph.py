'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import networkx as nx

from taco.lib.base import Exon, Strand
from taco.lib.transfrag import Transfrag

from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import choose_k, create_path_graph, get_kmers, \
    add_path, SOURCE, SINK, smooth_graph

from taco.test.base import read_single_locus


def test_choose_k_1():
    A = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 100), Exon(200, 300), Exon(400, 500)],
                  _id='A')
    SG = SpliceGraph.create([A])
    node_bounds = SG.node_bounds
    # case does not reach min frag length
    k = choose_k([A], node_bounds, min_path_length=400, kmin=1)
    assert k == 3
    # case reaches min frag length
    k = choose_k([A], node_bounds, min_path_length=300, kmin=1)
    assert k == 3
    k = choose_k([A], node_bounds, min_path_length=250, kmin=1)
    assert k == 3
    # case reaches frag length with nodes to spare
    k = choose_k([A], node_bounds, min_path_length=200, kmin=1)
    assert k == 2
    # case reaches frag length with nodes to spare
    k = choose_k([A], node_bounds, min_path_length=100, kmin=1)
    assert k == 1
    # enforce kmin
    k = choose_k([A], node_bounds, min_path_length=100, kmin=2)
    assert k == 2
    k = choose_k([A], node_bounds, min_path_length=100, kmin=3)
    assert k == 3
    k = choose_k([A], node_bounds, min_path_length=100, kmin=5)
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
    SG = SpliceGraph.create([A, B])
    node_bounds = SG.node_bounds
    k = choose_k([A, B], node_bounds, min_path_length=200, kmin=1)
    assert k == 11
    k = choose_k([A, B], node_bounds, min_path_length=100, kmin=1)
    assert k == 10
    k = choose_k([A, B], node_bounds, min_path_length=1, kmin=1)
    assert k == 1
    k = choose_k([A, B], node_bounds, min_path_length=1, kmin=3)
    assert k == 3

    B = Transfrag(chrom='chr1', strand=Strand.POS,
                  exons=[Exon(0, 10),
                         Exon(65, 75),
                         Exon(500, 600),
                         Exon(1000, 1100)],
                  _id='B')
    SG = SpliceGraph.create([A, B])
    node_bounds = SG.node_bounds
    k = choose_k([A, B], node_bounds, min_path_length=200)
    assert k == 12
    k = choose_k([A, B], node_bounds, min_path_length=100)
    assert k == 11
    k = choose_k([A, B], node_bounds, min_path_length=50)
    assert k == 6
    k = choose_k([A, B], node_bounds, min_path_length=1, kmin=1)
    assert k == 1
    k = choose_k([A, B], node_bounds, min_path_length=1, kmin=2)
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
    SG = SpliceGraph.create([A])
    node_bounds = SG.node_bounds
    k = choose_k([A], node_bounds, min_path_length=400, kmin=1)
    assert k == 1
    k = choose_k([A], node_bounds, min_path_length=400, kmin=2)
    assert k == 2
    k = choose_k([A], node_bounds, min_path_length=20000, kmin=2)
    assert k == 2
    k = choose_k([A], node_bounds, min_path_length=100000, kmin=2)
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
    G1, lost_paths = create_path_graph(SG, k)
    G2 = nx.DiGraph()
    for path in paths:
        kmers = list(get_kmers(path, k))
        add_path(G2, kmers, 1.0)
    assert nx.is_isomorphic(G1, G2)


def test_path_graph2():
    t_dict, locus = read_single_locus('change_point2.gtf')
    sgraph = SpliceGraph.create(t_dict.values())

    # trivial case without additional stops or starts
    k = choose_k(sgraph.transfrags, sgraph.node_bounds)
    assert k == 1
    K, lost_paths = create_path_graph(sgraph, k)
    kmer_id_map = dict((v, k) for k, v in K.graph['id_kmer_map'].iteritems())
    assert len(lost_paths) == 0
    n = kmer_id_map[(Exon(0, 100),)]
    assert K.node[n]['expr'] == 12.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0

    # add a stop site
    sgraph.stop_sites.add(50)
    sgraph.recreate_graph()
    K, lost_paths = create_path_graph(sgraph, k=2)
    assert len(lost_paths) == 0
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
    sgraph.recreate_graph()
    K, lost_paths = create_path_graph(sgraph, k=2)
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
