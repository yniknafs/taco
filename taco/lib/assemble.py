'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import logging

from gtf import GTF
from base import Strand
from transfrag import Transfrag
from locus import Locus
from splice_graph import SpliceGraph

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def get_splice_graphs(locus):
    for strand in Strand.POS, Strand.NEG, Strand.NA:
        transfrags = locus.get_transfrags(strand)
        if len(transfrags) > 0:
            sg = SpliceGraph.create(locus.get_transfrags(strand))
            yield sg

    #     # get connected components of graph which represent independent genes
    #     # unconnected components are considered different genes
    #     Gsubs = nx.weakly_connected_component_subgraphs(G)
    #     # add components as separate transcript graphs
    #     strand_graphs = []
    #     node_subgraph_map = {}
    #     for i,Gsub in enumerate(Gsubs):
    #         for n in Gsub:
    #             node_subgraph_map[n] = i
    #         tg = TranscriptGraph(chrom, strand, Gsub)
    #         tg.partial_paths = collections.defaultdict(lambda: 0.0)
    #         strand_graphs.append(tg)
    #     # populate transcript graphs with partial paths
    #     for t in transcript_list:
    #         # get original transcript nodes and subtract trimmed nodes
    #         # convert to collapsed nodes and bin according to subgraph
    #         # TODO: intronic transcripts may be split into multiple pieces,
    #         # should we allow this?
    #         subgraph_node_map = collections.defaultdict(lambda: set())
    #         for n in split_exons(t, G.graph['boundaries']):
    #             n = Exon(*n)
    #             if n in trim_nodes:
    #                 continue
    #             cn = node_chain_map[n]
    #             subgraph_id = node_subgraph_map[cn]
    #             subgraph_node_map[subgraph_id].add(cn)
    #         # add transcript node/score pairs to subgraphs
    #         for subgraph_id, subgraph_nodes in subgraph_node_map.iteritems():
    #             subgraph_nodes = sorted(subgraph_nodes,
    #                                     key=operator.attrgetter('start'),
    #                                     reverse=(strand == NEG_STRAND))
    #             tg = strand_graphs[subgraph_id]
    #             tg.partial_paths[tuple(subgraph_nodes)] += t.score
    #     transcript_graphs.extend(strand_graphs)
    # # convert
    # for tg in transcript_graphs:
    #     tg.partial_paths = tg.partial_paths.items()
    # return transcript_graphs


def assemble(gtf_file, output_dir,
             guided_strand=False,
             guided_ends=False,
             guided_assembly=False):
    # setup bedgraph output files
    file_prefix = os.path.join(output_dir, 'loci.unresolved')
    raw_bgfilehd = Locus.open_bedgraphs(file_prefix)
    file_prefix = os.path.join(output_dir, 'loci.resolved')
    resolved_bgfilehd = Locus.open_bedgraphs(file_prefix)

    # parse gtf file
    for interval, gtf_lines in GTF.parse_loci(open(gtf_file)):
        chrom, start, end = interval
        t_dict = Transfrag.parse_gtf(gtf_lines)

        locus = Locus.create(t_dict.values())
        logging.debug('Locus %s:%d-%d: '
                      '%d transfrags (+: %d, -: %d, .: %d)' %
                      (chrom, start, end, len(t_dict),
                       len(locus.strand_transfrags[Strand.POS]),
                       len(locus.strand_transfrags[Strand.NEG]),
                       len(locus.strand_transfrags[Strand.NA])))
        # write raw bedgraph files
        locus.write_bedgraph(raw_bgfilehd)
        # resolve unstranded transcripts
        num_resolved = locus.impute_unknown_strands()
        if num_resolved > 0:
            logging.debug('Locus %s:%d-%d: %d '
                          'resolved (+: %d, -: %d, .: %d)' %
                          (chrom, start, end, num_resolved,
                           len(locus.strand_transfrags[Strand.POS]),
                           len(locus.strand_transfrags[Strand.NEG]),
                           len(locus.strand_transfrags[Strand.NA])))
        # write bedgraph files after strand resolved
        locus.write_bedgraph(resolved_bgfilehd)

        # create splice graphs
        for sg in get_splice_graphs(locus):
            print 'splice graph', sg.chrom, Strand.to_gtf(sg.strand), sg.start, sg.end
            pass

    # close bedgraph files
    Locus.close_bedgraphs(raw_bgfilehd)
    Locus.close_bedgraphs(resolved_bgfilehd)
