#!/usr/bin/env python
'''
TACO: Transcriptome meta-assembly from RNA-Seq

@author: mkiyer
@author: yniknafs
'''
import sys
import logging

from taco.lib.run import Run

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

PROFILE = False
PROFILE_FILENAME = 'taco.prof'


def main():
    # instantiate from command line
    R = Run.create()
    R.log_args()
    logging.debug('Samples: %d' % (len(R.samples)))
    # merge gtf files
    msg = 'Aggregating GTF files'
    if R.status.aggregate:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.aggregate()
    # index loci
    msg = 'Indexing Loci'
    if R.status.index_loci:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.index_loci()
    # assemble
    msg = 'Assembling GTF files'
    if R.status.assemble:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.assemble()

    return 0


if __name__ == '__main__':
    if PROFILE:
        import cProfile
        sys.exit(cProfile.run('main()', filename=PROFILE_FILENAME))
    else:
        sys.exit(main())
