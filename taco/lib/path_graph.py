'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import array

from base import Exon, Strand
from graph import Graph
from optimize import maximize_bisect
from csuffixarray import SuffixArrayIndex

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def smooth_iteration(order, adjlist, exprs, smooth_exprs):
    smooths_tmp = array.array('f', [0.0 for x in xrange(len(smooth_exprs))])
    for i in order:
        smooth_expr = smooth_exprs[i]
        nbrs = adjlist[i]
        if len(nbrs) == 0:
            continue
        tot_nbr_expr = sum(exprs[j] for j in nbrs)
        if tot_nbr_expr == 0:
            # if all successors are zero apply smoothing evenly
            avg_expr = smooth_expr / len(nbrs)
            for j in nbrs:
                smooths_tmp[j] += avg_expr
                smooth_exprs[j] += avg_expr
        else:
            # apply smoothing proportionately
            for j in nbrs:
                frac_expr = (exprs[j] / tot_nbr_expr) * smooth_expr
                smooths_tmp[j] += frac_expr
                smooth_exprs[j] += frac_expr
    return smooths_tmp


def get_kmers(path, k):
    for i in xrange(0, len(path) - (k-1)):
        yield path[i:i+k]


class PathGraphFactory(object):
    def __init__(self, sgraph):
        self.chrom = sgraph.chrom
        self.start = sgraph.start
        self.end = sgraph.end
        self.strand = sgraph.strand
        self.transfrags = []
        self.paths = []
        self.exprs = []
        self.has_source = []
        self.has_sink = []
        start_nodes, stop_nodes = sgraph.get_start_stop_nodes()
        for t in sgraph.itertransfrags():
            path = sgraph.get_path(t)
            self.transfrags.append(t)
            self.paths.append(path)
            self.exprs.append(t.expr)
            self.has_source.append(path[0] in start_nodes)
            self.has_sink.append(path[-1] in stop_nodes)

    def longest_path_length(self):
        longest_path = 0
        for i in xrange(len(self.paths)):
            if len(self.paths[i]) > longest_path:
                longest_path = len(self.paths[i])
        return longest_path

    def create(self, k):
        short_transfrags = []
        K = PathGraph()
        for i in xrange(len(self.paths)):
            path = self.paths[i]
            full_length = self.has_source[i] and self.has_sink[i]
            short = len(path) < k
            if short and not full_length:
                short_transfrags.append(i)
                continue

            kmers = []
            if self.has_source[i]:
                kmers.append(K.SOURCE)
            if short:
                # specially add short full length paths because
                # they are not long enough to have kmers
                kmers.append(path)
            else:
                # get kmers
                kmers.extend(get_kmers(path, k))
            if self.has_sink[i]:
                kmers.append(K.SINK)
            K.add_path(kmers, self.exprs[i])

        num_lost_kmers = 0
        lost_kmer_expr = 0.0
        for i in K.get_unreachable_nodes():
            num_lost_kmers += 1
            lost_kmer_expr += K.exprs[i]
            K.remove_node_id(i)

        K.num_lost_kmers = num_lost_kmers
        K.lost_kmer_expr = lost_kmer_expr
        K.short_transfrags = short_transfrags
        K.valid = K.is_valid()
        return K

    def create_optimal(self, kmax=0, loss_threshold=0.10, stats_fh=None):
        '''
        create a graph where nodes are paths of length 'k'. the parameter
        'k' is chosen to maximizing the number of reachable k-mers in the
        path graph while tolerating at most 'loss_threshold' percent of
        expression.
        '''
        if len(self.paths) == 0:
            return None, 0

        # find upper bound to k
        user_kmax = kmax
        kmax = self.longest_path_length()
        if user_kmax > 0:
            # user can force a specific kmax (for debugging/testing purposes)
            kmax = min(user_kmax, kmax)
        id_str = ('%s:%d-%d[%s]' % (self.chrom, self.start, self.end,
                                    Strand.to_gtf(self.strand)))

        def compute_kmers(k):
            K = self.create(k)
            tot_expr = sum(K.exprs[i] for i in K.node_ids_iter())
            lost_expr = K.lost_kmer_expr
            lost_expr_frac = 0.0 if tot_expr == 0 else lost_expr / tot_expr
            logging.debug('%s k=%d kmax=%d t=%d kmers=%d short_transfrags=%d '
                          'lost_kmers=%d tot_expr=%.3f lost_expr=%.3f '
                          'lost_expr_frac=%.3f valid=%d' %
                          (id_str, k, kmax, len(self.paths), K.n,
                           len(K.short_transfrags), K.num_lost_kmers, tot_expr,
                           lost_expr, lost_expr_frac, int(K.valid)))
            if stats_fh:
                fields = [id_str, k, kmax, len(self.paths), K.n,
                          len(K.short_transfrags), K.num_lost_kmers, tot_expr,
                          lost_expr, lost_expr_frac, int(K.valid)]
                print >>stats_fh, '\t'.join(map(str, fields))
            if not K.valid:
                return -k
            #if lost_expr_frac > loss_threshold:
            #    return -k
            return len(K)

        k, num_kmers = maximize_bisect(compute_kmers, 1, kmax, 0)
        logging.debug('Creating path graph k=%d num_kmers=%d' % (k, num_kmers))
        K = self.create(k)
        logging.debug('Rescuing short transfrags kmers=%d' % len(K))
        self.rescue_short_transfrags(K, K.short_transfrags)
        return K, k

    def rescue_short_transfrags(self, K, indexes):
        # build suffix array index
        sai = K.get_suffix_array_index()
        # align short transfrags to index
        kmer_exprs = []
        lost_transfrags = []
        for i in indexes:
            # align to find matching kmers
            tot_expr = 0.0
            matching_kmers = []
            for kmer in sai.search(self.paths[i]):
                kmer_id = K.node_id_map[kmer]
                kmer_expr = K.exprs[kmer_id]
                tot_expr += kmer_expr
                matching_kmers.append((kmer_id, kmer_expr))
            # calculate expression for matching kmers
            matching_paths = []
            for kmer_id, kmer_expr in matching_kmers:
                if tot_expr == 0:
                    new_expr = self.exprs[i] / len(matching_kmers)
                else:
                    new_expr = self.exprs[i] * (kmer_expr / tot_expr)
                matching_paths.append((kmer_id, new_expr))
            if len(matching_paths) == 0:
                lost_transfrags.append(i)
                t = self.transfrags[i]
            kmer_exprs.extend(matching_paths)
        # add short transfrag kmers
        for kmer_id, expr in kmer_exprs:
            K.exprs[kmer_id] += expr
            K.smooth_rev[kmer_id] += expr
            K.smooth_fwd[kmer_id] += expr
        logging.debug('\tlost %d of %d short transfrags' %
                      (len(lost_transfrags), len(indexes)))


class PathGraph(Graph):
    def __init__(self):
        super(PathGraph, self).__init__()
        self.exprs = [0.0, 0.0]
        self.smooth_rev = [0.0, 0.0]
        self.smooth_fwd = [0.0, 0.0]

    def add_node(self, node):
        if node not in self.node_id_map:
            node_id = super(PathGraph, self).add_node(node)
            self.exprs.append(0.0)
            self.smooth_rev.append(0.0)
            self.smooth_fwd.append(0.0)
        else:
            node_id = self.node_id_map[node]
        return node_id

    def add_path(self, path, expr):
        node_ids = super(PathGraph, self).add_path(path)
        for node_id in node_ids:
            self.exprs[node_id] += expr
        # the first node should be "smoothed" in reverse direction
        self.smooth_rev[node_ids[0]] += expr
        # the last node should be "smoothed" in forward direction
        self.smooth_fwd[node_ids[-1]] += expr

    def apply_smoothing(self):
        order = self.topological_sort()
        # smooth in forward direction
        fwds = smooth_iteration(order, self.succs, self.exprs,
                                self.smooth_fwd)
        # smooth in reverse direction
        order.reverse()
        revs = smooth_iteration(order, self.preds, self.exprs,
                                self.smooth_rev)
        # apply densities to nodes
        for i in xrange(len(self.nodes)):
            self.exprs[i] += fwds[i] + revs[i]

    def get_suffix_array_index(self):
        # build suffix array index
        nodes = (self.nodes[i] for i in self.node_ids_iter())
        sai = SuffixArrayIndex(nodes)
        return sai


def reconstruct_path(kmer_path, pgraph, sgraph):
    # reconstruct nodes from kmers
    path = list(pgraph.nodes[kmer_path[1]])
    path.extend(pgraph.nodes[n][-1] for n in kmer_path[2:-1])
    # reverse negative stranded data so that all paths go from
    # small -> large genomic coords
    if sgraph.strand == Strand.NEG:
        path.reverse()
    # convert from integer node labels to genome (start, end) tuples
    path = [sgraph.get_node_interval(nid) for nid in path]
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
