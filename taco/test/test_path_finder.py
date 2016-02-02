'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

from taco.lib.base import Strand
from taco.lib.assemble import assemble_isoforms
from taco.lib.splice_graph import SpliceGraph

from taco.test.base import read_single_locus


def test_empty_graph_bug():
    t_dict, locus = read_single_locus('empty_graph_bug.gtf')
    transfrags = locus.get_transfrags(Strand.POS)
    sgraph = SpliceGraph.create(transfrags)
    isoforms = assemble_isoforms(sgraph,
                                 min_path_length=400,
                                 frac_isoform=0.01,
                                 max_isoforms=100)
    assert len(isoforms) == 0
