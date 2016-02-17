'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from cpython cimport array
from array import array
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
DEF MIN_SCORE = 1.0e-10


def topological_sort(object G):
    """From NetworkX source code
    https://github.com/networkx/networkx/blob/v1.10/networkx/algorithms/dag.py
    """
    cdef set explored
    cdef list order, fringe, new_nodes
    cdef int v, w, n

    order = []
    explored = set()

    for v in G.nodes_iter():     # process all vertices in G
        if v in explored:
            continue
        fringe = [v]   # nodes yet to look at
        while fringe:
            w = fringe[-1]  # depth first search
            if w in explored:  # already looked down this branch
                fringe.pop()
                continue
            # Check successors for cycles and for new nodes
            new_nodes = []
            for n in G[w]:
                if n not in explored:
                    new_nodes.append(n)
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored.add(w)
                order.append(w)
                fringe.pop()    # done considering this node
    order.reverse()
    return order


def find_path(int[:] nodes, float[:] exprs, list succ, int isource, int isink):
    cdef float[:] min_exprs
    cdef int[:] prevs
    cdef list ipath, path
    cdef int n, i, j, prev, x
    cdef float expr, min_expr, new_min_expr, new_expr

    cdef array.array int_array_template = array.array('i', [])
    cdef array.array float_array_template = array.array('f', [])

    # initialize data structures
    n = nodes.shape[0]
    min_exprs = array.clone(float_array_template, n, zero=False)
    prevs = array.clone(int_array_template, n, zero=False)
    for i in xrange(n):
        min_exprs[i] = MIN_SCORE
        prevs[i] = isink
    min_exprs[isource] = exprs[isource]

    # traverse nodes in topological sort order
    for i in xrange(n):
        min_expr = min_exprs[i]
        for j in succ[i]:
            new_min_expr = min_expr if min_expr < exprs[j] else exprs[j]
            if (prevs[j] == isink) or (new_min_expr > min_exprs[j]):
                min_exprs[j] = new_min_expr
                prevs[j] = i

    # traceback to get path
    expr = min_exprs[isink]
    prev = isink
    ipath = [isink]
    while True:
        prev = prevs[prev]
        ipath.append(prev)
        if prev == isource:
            break
    ipath.reverse()

    # subtract path
    for i in xrange(len(ipath)):
        x = ipath[i]
        new_expr = exprs[x] - expr
        exprs[x] = MIN_SCORE if MIN_SCORE >= new_expr else new_expr
        ipath[i] = nodes[x]
    return tuple(ipath), expr


def find_paths(object G, object expr_attr, float path_frac=0, int max_paths=0):
    cdef dict indexes
    cdef int[:] nodes
    cdef float[:] exprs
    cdef list succ
    cdef list results
    cdef int n, i, j, isource, isink, iterations
    cdef float expr, lowest_expr
    cdef tuple path

    # initialize data structures
    n = len(G)
    #nodes = array('i', nx.topological_sort(G))
    nodes = array('i', nx.topological_sort(G))
    indexes = {}
    exprs = array('f', nodes)
    succ = []
    for i in xrange(n):
        indexes[nodes[i]] = i
    for i in xrange(n):
        exprs[i] = G.node[nodes[i]][expr_attr]
        succ.append([indexes[x] for x in G.successors_iter(nodes[i])])
    isource = 0
    isink = n - 1

    # don't run if all nodes are zero
    if exprs[isource] < MIN_SCORE:
        return []

    # find highest scoring path
    path, expr = find_path(nodes, exprs, succ, isource, isink)
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
        path, expr = find_path(nodes, exprs, succ, isource, isink)
        if expr <= lowest_expr:
            break
        # store path
        results.append((path, expr))
        iterations += 1
    return results
