'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import logging

from gtf import GTF
from base import Strand
from transfrag import Transfrag
from locus import Locus

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def assemble(gtf_file,
             expr_h5_file,
             chrom_sizes_file,
             node_gtf_file,
             output_dir,
             guided_strand=False,
             guided_ends=False,
             guided_assembly=False):
    # setup bedgraph output files
    file_prefix = os.path.join(output_dir, 'loci.unresolved')
    raw_bgfilehd = Locus.open_bedgraphs(file_prefix)
    file_prefix = os.path.join(output_dir, 'loci.resolved')
    resolved_bgfilehd = Locus.open_bedgraphs(file_prefix)
    # setup expression hdf5
    expr_h5f = Locus.open_expression_hdf5(expr_h5_file, chrom_sizes_file)
    # setup node gtf file
    node_gtf_fileh = open(node_gtf_file, 'w')

    # parse gtf file
    for interval, gtf_lines in GTF.parse_loci(open(gtf_file)):
        chrom, start, end = interval
        t_dict = Transfrag.parse_gtf(gtf_lines)

        locus = Locus.create(t_dict.values())
        logging.debug('Locus %s:%d-%d: '
                      '%d transfrags (+: %d, -: %d, .: %d)' %
                      (chrom, start, end, len(t_dict),
                       len(locus.get_transfrags(Strand.POS)),
                       len(locus.get_transfrags(Strand.NEG)),
                       len(locus.get_transfrags(Strand.NA))))
        # write raw bedgraph files
        locus.write_bedgraph(raw_bgfilehd)

        # resolve unstranded transcripts
        num_resolved = locus.impute_unknown_strands()
        if num_resolved > 0:
            logging.debug('Locus %s:%d-%d: %d '
                          'resolved (+: %d, -: %d, .: %d)' %
                          (chrom, start, end, num_resolved,
                           len(locus.get_transfrags(Strand.POS)),
                           len(locus.get_transfrags(Strand.NEG)),
                           len(locus.get_transfrags(Strand.NA))))
        # write bedgraph files after strand resolved
        locus.write_bedgraph(resolved_bgfilehd)
        locus.write_expression_hdf5(expr_h5f)

        # convert to stranded locus objects
        for slocus in locus.create_stranded_loci():
            for f in slocus.get_node_gtf():
                print >>node_gtf_fileh, str(f)

    # close node file
    node_gtf_fileh.close()

    # close expression
    expr_h5f.close()

    # close bedgraph files
    Locus.close_bedgraphs(raw_bgfilehd)
    Locus.close_bedgraphs(resolved_bgfilehd)
