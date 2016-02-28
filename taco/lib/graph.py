'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def single_node_shortest_path_length(node, nbrs):
    '''
    Adapted from NetworkX 1.10 source code (networkx.github.io)

    node: source or sink
    nbrs: array of node neighbors (either succs or preds) indexed by node
    '''
    seen = {}
    level = 0
    nextlevel = set([node])
    while nextlevel:
        thislevel = nextlevel  # advance to next level
        nextlevel = set()      # and start a new list (fringe)
        for v in thislevel:
            if v not in seen:
                seen[v] = level  # set the level of vertex v
                nextlevel.update(nbrs[v])
        level += 1
    return seen  # return all path lengths as dictionary


class Graph(object):
    SOURCE = -1
    SINK = -2
    EMPTY = -3
    SOURCE_ID = 0
    SINK_ID = 1
    FIRST_ID = 2

    def __init__(self):
        self.node_id_map = {Graph.SOURCE: Graph.SOURCE_ID,
                            Graph.SINK: Graph.SINK_ID}
        self.nodes = [Graph.SOURCE, Graph.SINK]
        self.succs = [set(), set()]
        self.preds = [set(), set()]
        self.n = 0

    def __len__(self):
        return self.n

    def node_ids_iter(self, source=False, sink=False):
        if source:
            yield Graph.SOURCE_ID
        if sink:
            yield Graph.SINK_ID
        for i in xrange(Graph.FIRST_ID, len(self.nodes)):
            if self.nodes[i] == Graph.EMPTY:
                continue
            yield i

    def edges_iter(self):
        for i in xrange(len(self.nodes)):
            if self.nodes[i] == Graph.EMPTY:
                continue
            for j in self.succs[i]:
                yield i, j

    def has_node(self, node):
        return node in self.node_id_map

    def get_node(self, node):
        return self.node_id_map[node]

    def add_node(self, node):
        if node not in self.node_id_map:
            node_id = len(self.nodes)
            self.node_id_map[node] = node_id
            self.nodes.append(node)
            self.succs.append(set())
            self.preds.append(set())
            self.n += 1
        else:
            node_id = self.node_id_map[node]
        return node_id

    def add_path(self, nodes):
        if len(nodes) == 0:
            return []
        # add edges in path
        u = self.add_node(nodes[0])
        node_ids = [u]
        for i in xrange(1, len(nodes)):
            v = self.add_node(nodes[i])
            node_ids.append(v)
            self.succs[u].add(v)
            self.preds[v].add(u)
            u = v
        return node_ids

    def remove_node_id(self, node_id):
        if node_id == Graph.SOURCE_ID or node_id == Graph.SINK_ID:
            return
        # for each successor, remove from predecessors
        for succ in self.succs[node_id]:
            self.preds[succ].remove(node_id)
        # for each predecessor, remove from successors
        for pred in self.preds[node_id]:
            self.succs[pred].remove(node_id)
        # remove from node id map
        del self.node_id_map[self.nodes[node_id]]
        # mask
        self.nodes[node_id] = Graph.EMPTY
        self.n -= 1

    def get_unreachable_nodes(self):
        '''find unreachable kmers from source or sink'''
        allnodes = set(self.node_id_map.values())
        # unreachable from source
        a = set(single_node_shortest_path_length(Graph.SOURCE_ID, self.succs))
        a = allnodes - a
        # unreachable from sink
        b = set(single_node_shortest_path_length(Graph.SINK_ID, self.preds))
        b = allnodes - b
        return a | b

    def remove_unreachable_nodes(self):
        '''
        mask nodes that are unreachable from the source or sink, these occur
        due to fragmentation when k > 1
        '''
        unreachable = self.get_unreachable_nodes()
        for i in unreachable:
            self.remove_node_id(i)
        return len(unreachable)

    def is_valid(self):
        '''
        Adapted from NetworkX 1.10 bidirectional_shortest_path
        '''
        if self.nodes[Graph.SOURCE_ID] == Graph.EMPTY:
            return False
        if self.nodes[Graph.SINK_ID] == Graph.EMPTY:
            return False
        # predecesssor and successors in search
        pred = {Graph.SOURCE_ID: None}
        succ = {Graph.SINK_ID: None}
        # initialize fringes, start with forward
        forward_fringe = [Graph.SOURCE_ID]
        reverse_fringe = [Graph.SINK_ID]

        while forward_fringe and reverse_fringe:
            if len(forward_fringe) <= len(reverse_fringe):
                this_level = forward_fringe
                forward_fringe = []
                for v in this_level:
                    for w in self.succs[v]:
                        if w not in pred:
                            forward_fringe.append(w)
                            pred[w] = v
                        if w in succ:
                            return True
                            # return pred, succ, w  # found path
            else:
                this_level = reverse_fringe
                reverse_fringe = []
                for v in this_level:
                    for w in self.preds[v]:
                        if w not in succ:
                            succ[w] = v
                            reverse_fringe.append(w)
                        if w in pred:
                            return True
                            # return pred, succ, w  # found path
        return False

    def topological_sort_kahn(self):
        '''
        Adapted from NetworkX source code (networkx.github.io)
        networkx.algorithms.dag.topological_sort
        '''
        indegree_map = {}
        zero_indegree = []
        for i in xrange(len(self.preds)):
            if self.nodes[i] == Graph.EMPTY:
                continue
            indegree = len(self.preds[i])
            if indegree > 0:
                indegree_map[i] = indegree
            else:
                zero_indegree.append(i)

        nodes = []
        while zero_indegree:
            i = zero_indegree.pop()
            for child in self.succs[i]:
                indegree_map[child] -= 1
                if indegree_map[child] == 0:
                    zero_indegree.append(child)
                    del indegree_map[child]
            nodes.append(i)

        if indegree_map:
            return None
        return nodes

    def topological_sort_dfs(self):
        """
        Adapted from NetworkX source code (networkx.github.io)
        networkx.algorithms.dag.topological_sort
        """
        # nonrecursive version
        order = []
        explored = set()
        for i in xrange(len(self.nodes)):
            if self.nodes[i] == Graph.EMPTY:
                continue
            if i in explored:
                continue
            fringe = [i]   # nodes yet to look at
            while fringe:
                j = fringe[-1]  # depth first search
                if j in explored:  # already looked down this branch
                    fringe.pop()
                    continue
                # Check successors for cycles and for new nodes
                new_nodes = [n for n in self.succs[j] if n not in explored]
                if new_nodes:   # Add new_nodes to fringe
                    fringe.extend(new_nodes)
                else:           # No new nodes so w is fully explored
                    explored.add(j)
                    order.append(j)
                    fringe.pop()    # done considering this node
        order.reverse()
        return order

    def topological_sort(self):
        return self.topological_sort_dfs()

    def is_topological_sort(self, order):
        order_indexes = dict((x, i) for i, x in enumerate(order))
        a = set(self.node_id_map.values()).symmetric_difference(order_indexes)
        if len(a) != 0:
            return False

        for u, v in self.edges_iter():
            ui = order_indexes[u]
            vi = order_indexes[v]
            if ui > vi:
                return False
        return True
