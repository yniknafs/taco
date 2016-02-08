'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from operator import itemgetter

from taco.lib.assemble import Cluster


def test_cluster():
    a1 = ((1, 2, 3, 4), 1000)
    a2 = ((4, 5, 6), 1000)
    b1 = ((7, 8, 9), 100)
    rt1 = ((6, 7), 1)

    clusters, filtered = Cluster.build([a1, a2, b1, rt1], min_frac=0.0)
    assert len(clusters) == 1
    assert len(filtered) == 0

    clusters, filtered = Cluster.build([a1, a2, b1, rt1], min_frac=0.01)
    assert len(clusters) == 2
    assert len(filtered) == 1
    # print 'clusters', len(clusters)
    # for c in clusters:
    #     print '_id', c._id, 'expr', c.expr, 'nodes', c.nodes, 'paths', c.paths
    # print 'filtered', len(filtered), filtered
