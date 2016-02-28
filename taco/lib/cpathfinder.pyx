'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from cpython cimport array
import array

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


def find_path(int[:] order, float[:] exprs, list succs,
              int source, int sink):
    cdef float[:] min_exprs
    cdef float[:] sum_exprs
    cdef int[:] lengths
    cdef int[:] prevs
    cdef float min_expr, sum_expr, new_min_expr, new_sum_expr
    cdef float new_avg_expr, cur_avg_expr
    cdef int n, i, j, length, new_length, prev
    cdef list path
    cdef array.array int_array_template = array.array('i', [])
    cdef array.array float_array_template = array.array('f', [])

    # initialize data structures
    n = exprs.shape[0]
    min_exprs = array.clone(float_array_template, n, zero=False)
    sum_exprs = array.clone(float_array_template, n, zero=False)
    lengths = array.clone(int_array_template, n, zero=False)
    prevs = array.clone(int_array_template, n, zero=False)
    for i in xrange(n):
        min_exprs[i] = MIN_SCORE
        sum_exprs[i] = MIN_SCORE
        lengths[i] = 1
        prevs[i] = sink
    min_exprs[source] = exprs[source]
    sum_exprs[source] = exprs[source]

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
#            if ((prevs[j] == sink) or (new_min_expr > min_exprs[j]) or
#                (new_min_expr == min_exprs[j] and
#                 new_avg_expr > cur_avg_expr)):
            if ((prevs[j] == sink) or (new_min_expr > min_exprs[j])):
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
        j = path[i]
        new_expr = exprs[j] - expr
        exprs[j] = MIN_SCORE if MIN_SCORE >= new_expr else new_expr
    return tuple(path), expr


def find_paths(object G, float path_frac=0, int max_paths=0):
    cdef int[:] order
    cdef float[:] exprs
    cdef list results
    cdef int n, i, j, source, sink, iterations
    cdef float expr, lowest_expr
    cdef tuple path

    # don't run if all nodes are zero
    if G.exprs[G.SOURCE_ID] < MIN_SCORE:
        return []

    # initialize data structures
    order = array.array('i', G.topological_sort())
    exprs = array.array('f', G.exprs)
    source = G.SOURCE_ID
    sink = G.SINK_ID

    # find highest scoring path
    path, expr = find_path(order, exprs, G.succs, source, sink)
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
        path, expr = find_path(order, exprs, G.succs, source, sink)
        if expr <= lowest_expr:
            break
        # store path
        results.append((path, expr))
        iterations += 1
    return results
