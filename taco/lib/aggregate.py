'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer
'''
import os
import logging
import collections
import pysam

from base import Sample, TacoError
from gtf import GTF, GTFError, sort_gtf
from stats import scoreatpercentile
from caggregate import caggregate

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "mkiyer@umich.edu"
__status__ = "Development"


def aggregate(samples, ref_gtf_file, gtf_expr_attr, tmp_dir,
              output_gtf_file, stats_file):
    '''
    Aggregate/merge individual sample GTF files
    '''
    # setup output files
    tmp_file = os.path.join(tmp_dir, 'transcripts.unsorted.gtf')

    # aggregate ref gtf
    if ref_gtf_file is not None:
        logging.debug('Reference: %s' % ref_gtf_file)
        caggregate(ref_gtf_file, str(Sample.REF_ID), gtf_expr_attr,
                   tmp_file, stats_file, str(True))

    # aggregate sample gtfs
    for sample in samples:
        logging.debug('Sample: %s %s' % (sample._id, sample.gtf_file))
        caggregate(sample.gtf_file, str(sample._id), gtf_expr_attr,
                   tmp_file, stats_file, str(False))

    # sort merged gtf
    logging.info("Sorting GTF")
    retcode = sort_gtf(tmp_file, output_gtf_file, tmp_dir=tmp_dir)
    if retcode != 0:
        logging.error("Error sorting GTF")
        if os.path.exists(output_gtf_file):
            os.remove(output_gtf_file)
        raise TacoError('Error sorting GTF')
    os.remove(tmp_file)
