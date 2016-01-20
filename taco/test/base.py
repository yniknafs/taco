'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os

from taco.lib.gtf import GTF
from taco.lib.transfrag import Transfrag
from taco.lib.assemble import Locus

INPUT_FILE_DIR = "input_files"


def get_gtf_path(filename):
    return os.path.join(os.path.dirname(__file__), INPUT_FILE_DIR, filename)


def read_gtf(filename):
    return list(GTF.parse_loci(open(get_gtf_path(filename))))


def read_single_locus(filename):
    loci = read_gtf(filename)
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    locus = Locus.create(t_dict.values())
    return locus
