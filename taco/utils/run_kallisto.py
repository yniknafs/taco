#!/usr/bin/env python

'''
@author yniknafs
'''

import os
import sys
import argparse
import logging
import collections
import subprocess
import glob



def main():
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin', dest='bin',
            default = '/mctp/projects/taco/sw/kallisto_linux-v0.42.4/kallisto')
    parser.add_argument('-i', dest='index',
                        default='/mctp/projects/taco/ref/kallisto_index_gencode_v22')
    parser.add_argument('--fq1', dest='fq1')
    parser.add_argument('--fq2', dest='fq2')
    parser.add_argument('-p', dest='procs', default='1')
    parser.add_argument('out')
    args = parser.parse_args()


    args = [
    args.bin,
    'quant',
    '-i', args.index,
    '-o', args.out,
    '-t', args.procs,
    args.fq1,
    args.fq2,
    ]

    cmd = [
    '/bin/bash', '-c', ' '.join(args)
    ]

    subprocess.call(cmd)

    return 0

if __name__ == '__main__':
    sys.exit(main())
