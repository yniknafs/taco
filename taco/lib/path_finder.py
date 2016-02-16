'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import networkx as nx

from path_graph import KMER_EXPR

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


# constant minimum path score
MIN_SCORE = 1.0e-10


def _subtract_path(ipath, expr, exprs):
    for i in ipath:
        new_expr = exprs[i] - expr
        exprs[i] = MIN_SCORE if MIN_SCORE >= new_expr else new_expr


def _find_path(nodes, exprs, succ, isource, isink):
    min_exprs = [MIN_SCORE for i in xrange(len(exprs))]
    min_exprs[isource] = exprs[isource]
    prevs = [None for i in xrange(len(exprs))]
    for i in xrange(len(exprs)):
        min_expr = min_exprs[i]
        for j in succ[i]:
            new_min_expr = min_expr if min_expr < exprs[j] else exprs[j]
            if (prevs[j] is None) or (new_min_expr > min_exprs[j]):
                min_exprs[j] = new_min_expr
                prevs[j] = i
    # traceback to get path
    expr = min_exprs[isink]
    ipath = [isink]
    prev = prevs[isink]
    while prev is not None:
        ipath.append(prev)
        prev = prevs[prev]
    ipath.reverse()
    _subtract_path(ipath, expr, exprs)
    path = tuple(nodes[i] for i in ipath)
    return path, expr


def find_paths(G, source, sink, path_frac=0, max_paths=0, relative_frac=False):
    # initialize path finding
    nodes = nx.topological_sort(G)
    indexes = dict((n, i) for i, n in enumerate(nodes))
    isource, isink = 0, len(nodes)-1
    exprs = []
    succ = []
    for n in nodes:
        exprs.append(G.node[n][KMER_EXPR])
        succ.append([indexes[x] for x in G.successors_iter(n)])

    # don't run if all nodes are zero
    if exprs[isource] < MIN_SCORE:
        return []

    # find highest scoring path
    path, expr = _find_path(nodes, exprs, succ, isource, isink)
    results = [(path, expr)]

    # define threshold score to stop producing suboptimal paths
    if relative_frac:
        lowest_expr = expr * path_frac
    else:
        lowest_expr = exprs[isource] * path_frac
    lowest_expr = max(MIN_SCORE, lowest_expr)

    # iterate to find suboptimal paths
    iterations = 1
    while True:
        if max_paths > 0 and iterations >= max_paths:
            break
        # find path
        path, expr = _find_path(nodes, exprs, succ, isource, isink)
        if expr <= lowest_expr:
            break
        # store path
        results.append((path, expr))
        iterations += 1
    logging.debug("\tpath finding iterations=%d" % iterations)
    # return (path,score) tuples sorted from high -> low score
    return results
