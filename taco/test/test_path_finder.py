'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

from taco.lib.base import Strand
from taco.lib.assemble import assemble_isoforms, Config
from taco.lib.splice_graph import SpliceGraph

from taco.test.base import read_single_locus


def test_empty_graph_bug():
    t_dict, locus = read_single_locus('empty_graph_bug.gtf')
    transfrags = locus.get_transfrags(Strand.POS)
    sgraph = SpliceGraph.create(transfrags)
    isoforms = assemble_isoforms(sgraph, Config.defaults())
    assert len(isoforms) == 0
