'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import logging

from taco.lib.gtf import GTF
from taco.lib.transfrag import Transfrag
from taco.lib.assemble import find_boundaries, split_transfrag, Locus

INPUT_FILE_DIR = "input_files"


def get_gtf_path(filename):
    return os.path.join(os.path.dirname(__file__), INPUT_FILE_DIR, filename)


def read_gtf(filename):
    return list(GTF.parse_loci(open(get_gtf_path(filename))))


def test_parse_loci():
    loci = read_gtf('parse_loci.gtf')
    assert len(loci) == 3
    assert loci[0][0] == ('chrTest1', 10, 50)
    assert loci[1][0] == ('chrTest1', 50, 200)
    assert loci[2][0] == ('chrTest2', 100, 200)


def test_find_splice_sites():
    loci = read_gtf('splice_sites.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 10, 525)
    t_dict = Transfrag.parse_gtf(gtf_lines)
    splice_sites = tuple(find_boundaries(t_dict.values()))
    assert splice_sites == (10, 100, 200, 250, 300, 400, 525)


def test_find_transfrag_nodes():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    boundaries = find_boundaries(t_dict.values())
    # check nodes
    t = t_dict['A']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (200, 250), (250, 300), (400, 525))
    t = t_dict['B']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((10, 100), (250, 300), (400, 525))
    t = t_dict['C']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((100, 200), (200, 250), (250, 300), (400, 525))
    t = t_dict['D']
    nodes = tuple(split_transfrag(t, boundaries))
    assert nodes == ((300, 400), (400, 525))


def test_create_locus():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    locus = Locus.create(t_dict.values())
    # logging.error()
