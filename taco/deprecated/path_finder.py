'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import logging
import networkx as nx


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


def find_paths(G, expr_attr, path_frac=0, max_paths=0):
    # initialize path finding
    nodes = nx.topological_sort(G)
    indexes = dict((n, i) for i, n in enumerate(nodes))
    isource, isink = 0, len(nodes)-1
    exprs = []
    succ = []
    for n in nodes:
        exprs.append(G.node[n][expr_attr])
        succ.append([indexes[x] for x in G.successors_iter(n)])

    # don't run if all nodes are zero
    if exprs[isource] < MIN_SCORE:
        return []
    # find highest scoring path
    path, expr = _find_path(nodes, exprs, succ, isource, isink)
    results = [(path, expr)]
    # define threshold score to stop producing suboptimal paths
    lowest_expr = expr * path_frac
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
    return results


def find_path2(order, exprs, succs, source, sink):
    # initialize data structures
    min_exprs = []
    sum_exprs = []
    lengths = []
    prevs = []
    for i in xrange(len(exprs)):
        min_exprs.append(MIN_SCORE)
        sum_exprs.append(MIN_SCORE)
        lengths.append(1)
        prevs.append(sink)
    min_exprs[source] = exprs[source]
    sum_exprs[source] = exprs[source]
    lengths[source] = 1

    # traverse nodes in topological sort order
    for i in order:
        min_expr = min_exprs[i]
        sum_expr = sum_exprs[i]
        length = lengths[i]

        for j in succs[i]:
            new_min_expr = min_expr if min_expr < exprs[j] else exprs[j]
            new_sum_expr = sum_expr + exprs[j]
            new_length = length + 1
            new_avg_expr = new_sum_expr / new_length
            cur_avg_expr = sum_exprs[j] / lengths[j]

            # update node if 1) not yet visited, 2) can reach this node
            # with a higher min expr, or 3) can reach with an equal
            # min expr but a higher overall average expr
            update = ((prevs[j] == sink) or
                      (new_min_expr > min_exprs[j]) or
                      (new_min_expr == min_exprs[j] and
                       new_avg_expr > cur_avg_expr))
            if update:
                min_exprs[j] = new_min_expr
                sum_exprs[j] = new_sum_expr
                lengths[j] = new_length
                prevs[j] = i

    # traceback to get path
    expr = min_exprs[sink]
    prev = sink
    path = [sink]
    while True:
        prev = prevs[prev]
        path.append(prev)
        if prev == source:
            break
    path.reverse()

    # subtract path
    for i in xrange(len(path)):
        x = path[i]
        new_expr = exprs[x] - expr
        exprs[x] = MIN_SCORE if MIN_SCORE >= new_expr else new_expr
    return tuple(path), expr


def find_paths2(G, path_frac=0, max_paths=0):
    # initialize data structures
    assert G.is_topological_sort(G.topological_sort())
    order = G.topological_sort()
    exprs = list(G.exprs)

    # don't run if all nodes are zero
    if G.exprs[G.SOURCE_ID] < MIN_SCORE:
        return []
    # find highest scoring path
    path, expr = find_path2(order, exprs, G.succs, G.SOURCE_ID, G.SINK_ID)
    results = [(path, expr)]

    # define threshold score to stop producing paths
    lowest_expr = expr * path_frac
    if MIN_SCORE > lowest_expr:
        lowest_expr = MIN_SCORE

    # iterate to find paths
    iterations = 1
    while True:
        if max_paths > 0 and iterations >= max_paths:
            break
        # find path
        path, expr = find_path2(order, exprs, G.succs, G.SOURCE_ID, G.SINK_ID)
        if expr <= lowest_expr:
            break
        # store path
        results.append((path, expr))
        iterations += 1
    return results
