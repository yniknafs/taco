#!/usr/bin/env python
'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import sys
import argparse
import logging
import subprocess

from taco.lib.run import Results

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def prepare_igvtools(igvtools, infile, ext='gtf'):
    outfile = os.path.splitext(infile)[0] + '.srt.' + ext
    logging.info('igvtools sort %s -> %s' % (infile, outfile))
    subprocess.call([igvtools, 'sort', infile, outfile])
    logging.info('igvtools index %s' % (outfile))
    subprocess.call([igvtools, 'index', outfile])


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedGraphToBigWig', dest='bedGraphToBigWig',
                        help='path to "bedGraphToBigWig" binary executable')
    parser.add_argument('--igvtools', dest='igvtools',
                        help='path to "igvtools" binary executable')
    parser.add_argument('--chrom-sizes-file', dest='chrom_sizes_file',
                        help='path to chromosome sizes file correponding to '
                        'genome used to align reads')
    parser.add_argument('taco_dir',
                        help='path to taco run output directory')
    args = parser.parse_args()

    bedGraphToBigWig = None
    if args.bedGraphToBigWig and is_exe(args.bedGraphToBigWig):
        bedGraphToBigWig = args.bedGraphToBigWig
    else:
        bedGraphToBigWig = which('bedGraphToBigWig')
    if bedGraphToBigWig is None:
        parser.error('bedGraphToBigWig executable not found')

    igvtools = None
    if args.igvtools and is_exe(args.igvtools):
        igvtools = args.igvtools
    else:
        igvtools = which('igvtools')
    if igvtools is None:
        parser.error('igvtools executable not found')

    chrom_sizes_file = None
    if args.chrom_sizes_file and os.path.exists(args.chrom_sizes_file):
        chrom_sizes_file = args.chrom_sizes_file
    else:
        parser.error('chromosome sizes file not found')

    output_dir = args.taco_dir
    if not os.path.exists(output_dir):
        parser.error('TACO output directory %s not found' % output_dir)
    if not os.path.isdir(output_dir):
        parser.error('%s is not a directory' % output_dir)
    r = Results(output_dir)

    logging.info('bedGraphToBigWig: %s' % bedGraphToBigWig)
    logging.info('igvtools: %s' % igvtools)
    logging.info('chrom_sizes_file: %s' % chrom_sizes_file)
    logging.info('output directory %s' % output_dir)

    for strand, infile in enumerate(r.bedgraph_files):
        outfile = os.path.splitext(infile)[0] + '.bw'
        logging.info('Converting bedGraph to bigWig %s -> %s' %
                     (infile, outfile))
        subprocess.call([bedGraphToBigWig, infile, chrom_sizes_file,
                         outfile])

    prepare_igvtools(igvtools, r.splice_graph_gtf_file, ext='gtf')
    prepare_igvtools(igvtools, r.assembly_gtf_file, ext='gtf')
    prepare_igvtools(igvtools, r.splice_bed_file, ext='bed')


if __name__ == '__main__':
    sys.exit(main())
