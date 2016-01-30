'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import collections
import networkx as nx

from base import Exon, Strand
from splice_graph import split_transfrag

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


# graph attributes
SOURCE = -1
SINK = -2
NODE_EXPR = 'expr'
SMOOTH_FWD = 'smfwd'
SMOOTH_REV = 'smrev'
SMOOTH_TMP = 'smtmp'


def smooth_iteration(G, expr_attr, smooth_attr):
    nodes = nx.topological_sort(G)
    for u in nodes:
        ud = G.node[u]
        smooth_expr = ud[smooth_attr]
        succ = G.successors(u)
        if len(succ) == 0:
            continue
        total_nbr_expr = sum(G.node[v][expr_attr] for v in succ)
        if total_nbr_expr == 0:
            # if all successors have zero score apply smoothing evenly
            avg_expr = smooth_expr / len(succ)
            for v in succ:
                vd = G.node[v]
                vd[SMOOTH_TMP] += avg_expr
                vd[smooth_attr] += avg_expr
        else:
            # apply smoothing proportionately
            for v in succ:
                vd = G.node[v]
                frac = vd[expr_attr]/float(total_nbr_expr)
                adj_expr = frac * smooth_expr
                vd[SMOOTH_TMP] += adj_expr
                vd[smooth_attr] += adj_expr


def smooth_graph(G, expr_attr=NODE_EXPR):
    # smooth in forward direction
    smooth_iteration(G, expr_attr, SMOOTH_FWD)
    # smooth in reverse direction
    G.reverse(copy=False)
    smooth_iteration(G, expr_attr, SMOOTH_REV)
    G.reverse(copy=False)
    # apply densities to nodes
    for n, d in G.nodes_iter(data=True):
        d[expr_attr] += d[SMOOTH_TMP]


def choose_k(transfrags, node_bounds, min_path_length=400, kmin=2):
    assert kmin > 0
    k = 1
    for t in transfrags:
        node_lengths = [(n[1] - n[0]) for n in split_transfrag(t, node_bounds)]
        # initialize
        path_length = 0
        path_k = 0
        i = 0
        j = 0
        while i < len(node_lengths):
            while j < len(node_lengths):
                if (path_length >= min_path_length) and ((j - i) >= kmin):
                    break
                path_length += node_lengths[j]
                j += 1
            path_k = max(path_k, (j - i))
            if j == len(node_lengths):
                break
            path_length -= node_lengths[i]
            i += 1
        k = max(k, path_k)
    return k


def init_node_attrs():
    return {NODE_EXPR: 0.0, SMOOTH_FWD: 0.0,
            SMOOTH_REV: 0.0, SMOOTH_TMP: 0.0}


def add_path(K, path, expr):
    # add first kmer
    from_id = path[0]
    if from_id not in K:
        K.add_node(from_id, attr_dict=init_node_attrs())
    kmerattrs = K.node[from_id]
    kmerattrs[NODE_EXPR] += expr
    # the first kmer should be "smoothed" in reverse direction
    kmerattrs[SMOOTH_REV] += expr
    for to_id in path[1:]:
        if to_id not in K:
            K.add_node(to_id, attr_dict=init_node_attrs())
        kmerattrs = K.node[to_id]
        kmerattrs[NODE_EXPR] += expr
        # connect kmers
        K.add_edge(from_id, to_id)
        # update from_kmer to continue loop
        from_id = to_id
    # the last kmer should be "smoothed" in forward direction
    kmerattrs[SMOOTH_FWD] += expr


def hash_kmers(id_kmer_map, k, ksmall):
    kmer_hash = collections.defaultdict(lambda: set())
    for kmer_id, kmer in id_kmer_map.iteritems():
        for i in xrange(k - (ksmall-1)):
            kmer_hash[kmer[i:i+ksmall]].add(kmer_id)
    return kmer_hash


def find_short_path_kmers(kmer_hash, K, path, expr):
    """
    find kmers where 'path' is a subset and partition 'expr'
    of path proportionally among all matching kmers

    generator function yields (kmer_id, expr) tuples
    """
    if path not in kmer_hash:
        return
    matching_kmers = []
    total_expr = 0.0
    for kmer_id in kmer_hash[path]:
        # compute total expr at matching kmers
        kmer_expr = K.node[kmer_id][NODE_EXPR]
        total_expr += kmer_expr
        matching_kmers.append((kmer_id, kmer_expr))
    # now calculate fractional densities for matching kmers
    for kmer_id, kmer_expr in matching_kmers:
        new_expr = expr * (kmer_expr / float(total_expr))
        yield ([kmer_id], new_expr)


def connect_dangling_ends(K):
    """
    fragmented transcripts can manifest as 0-degree dangling ends
    of the overlap graph. this function connects these ends to the
    'source' and/or 'sink' nodes
    """
    source = K.graph['source']
    sink = K.graph['sink']
    # connect all nodes with degree zero to the source/sink nodes
    # to account for fragmentation in the kmer graph when k > 2
    edges_to_add = []
    for n in K.nodes_iter():
        if (n == source) or (n == sink):
            continue
        if K.in_degree(n) == 0:
            edges_to_add.append((source, n))
        if K.out_degree(n) == 0:
            edges_to_add.append((n, sink))
    # connect kmers
    for u, v in edges_to_add:
        K.add_edge(u, v)


def get_kmers(path, k):
    for i in xrange(0, len(path) - (k-1)):
        yield path[i:i+k]


def get_path(sgraph, t):
    nodes = [Exon(*n) for n in split_transfrag(t, sgraph.node_bounds)]
    if sgraph.strand == Strand.NEG:
        nodes.reverse()
    return tuple(nodes)


def create_path_graph(sgraph, k):
    '''create kmer graph from partial paths'''
    # initialize path graph
    K = nx.DiGraph()
    K.graph['source'] = SOURCE
    K.graph['sink'] = SINK
    K.add_node(SOURCE, attr_dict=init_node_attrs())
    K.add_node(SINK, attr_dict=init_node_attrs())
    # find all beginning/end nodes in splice graph
    start_nodes, stop_nodes = sgraph.get_start_stop_nodes()
    # convert paths to k-mers and create a k-mer to integer node map
    kmer_id_map = {}
    id_kmer_map = {}
    kmer_paths = []
    current_id = 0
    short_partial_path_dict = collections.defaultdict(lambda: [])
    for t in sgraph.get_transfrags():
        # get nodes
        path = get_path(sgraph, t)
        # check for start and stop nodes
        is_start = (path[0] in start_nodes)
        is_end = (path[-1] in stop_nodes)
        full_length = is_start and is_end
        if (len(path) < k) and (not full_length):
            # save fragmented short paths
            short_partial_path_dict[len(path)].append((path, t.expr))
            continue
        # convert to path of kmers
        kmer_path = []
        if is_start:
            kmer_path.append(SOURCE)
        if len(path) < k:
            # specially add short full length paths because
            # they are not long enough to have kmers
            kmers = [path]
        else:
            kmers = get_kmers(path, k)
        # convert to path of kmers
        for kmer in kmers:
            if kmer not in kmer_id_map:
                kmer_id = current_id
                kmer_id_map[kmer] = current_id
                id_kmer_map[current_id] = kmer
                current_id += 1
            else:
                kmer_id = kmer_id_map[kmer]
            kmer_path.append(kmer_id)
        if is_end:
            kmer_path.append(SINK)
        kmer_paths.append((kmer_path, t.expr))
    # store mapping from kmer_id to subpath tuple
    K.graph['id_kmer_map'] = id_kmer_map
    for path, expr in kmer_paths:
        add_path(K, path, expr)
    # try to add short paths to graph if they are exact subpaths of
    # existing kmers
    kmer_paths = []
    lost_paths = []
    for ksmall, short_partial_paths in short_partial_path_dict.iteritems():
        kmer_hash = hash_kmers(id_kmer_map, k, ksmall)
        for path, expr in short_partial_paths:
            matching_paths = \
                list(find_short_path_kmers(kmer_hash, K, path, expr))
            if len(matching_paths) == 0:
                lost_paths.append((path, expr))
            kmer_paths.extend(matching_paths)
    # add new paths
    for path, expr in kmer_paths:
        add_path(K, path, expr)
    # connect all kmer nodes with degree zero to the source/sink node
    # to account for fragmentation in the kmer graph when k > 2
    connect_dangling_ends(K)
    # cleanup
    del kmer_id_map
    return K, lost_paths


def reconstruct_path(kmer_path, id_kmer_map, strand):
    # reconstruct path from kmer ids
    path = list(id_kmer_map[kmer_path[1]])
    path.extend(id_kmer_map[n][-1] for n in kmer_path[2:-1])
    # reverse negative stranded data so that all paths go from
    # small -> large genomic coords
    if strand == Strand.NEG:
        path.reverse()
    # collapse contiguous nodes along path
    newpath = []
    chain = [path[0]]
    for v in path[1:]:
        if chain[-1].end != v.start:
            # update path with merge chain node
            newpath.append(Exon(chain[0].start,
                                chain[-1].end))
            # reset chain
            chain = []
        chain.append(v)
    # add last chain
    newpath.append(Exon(chain[0].start, chain[-1].end))
    return newpath
