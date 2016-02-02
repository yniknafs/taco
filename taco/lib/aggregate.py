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

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "mkiyer@umich.edu"
__status__ = "Development"


def _aggregate_gtf(gtf_file, sample_id, gtf_expr_attr, output_fh, stats_fh,
                   is_ref=False):
    def _init_t_dict():
        return {'_id': None, 'num_exons': 0, 'length': 0}

    t_dict = collections.defaultdict(_init_t_dict)
    cur_t_id = 1
    exprs = []
    for f in pysam.tabix_iterator(open(gtf_file), pysam.asGTF()):
        if f.feature == 'transcript':
            t_id = f.transcript_id
            if t_id in t_dict:
                m = 'GTF "%s" transcript_id "%s" not unique' % (gtf_file, t_id)
                raise GTFError(m)
            t_item = t_dict[t_id]
            # rename transcript id
            new_t_id = "%s.T%d" % (sample_id, cur_t_id)
            cur_t_id += 1
            t_item['_id'] = new_t_id
            if is_ref:
                expr = 0.0
            else:
                expr = float(f[gtf_expr_attr])
            exprs.append(expr)
            # prepare attributes
            attrs = {GTF.Attr.TRANSCRIPT_ID: new_t_id,
                     GTF.Attr.SAMPLE_ID: sample_id,
                     GTF.Attr.REF: str(int(is_ref)),
                     GTF.Attr.EXPR: str(expr)}
            # save attributes
            f.fromDict(attrs)
            print >>output_fh, str(f)
        elif f.feature == 'exon':
            t_id = f.transcript_id
            t_item = t_dict[t_id]
            # update statistics
            t_item['num_exons'] += 1
            t_item['length'] += (f.end - f.start)
            # replace transcript id
            f.fromDict({GTF.Attr.TRANSCRIPT_ID: t_item['_id']})
            print >>output_fh, str(f)

    # process statistics
    num_exons = []
    lengths = []
    for t_item in t_dict.itervalues():
        lengths.append(t_item['length'])
        num_exons.append(t_item['num_exons'])

    # compute and write stats
    quantiles = range(0, 101)
    expr_qs = (scoreatpercentile(exprs, q) for q in quantiles)
    expr_qs = ','.join(map(str, expr_qs))
    length_qs = (int(round(scoreatpercentile(lengths, q)))
                 for q in quantiles)
    length_qs = ','.join(map(str, length_qs))
    num_exon_qs = (int(round(scoreatpercentile(num_exons, q)))
                   for q in quantiles)
    num_exon_qs = ','.join(map(str, num_exon_qs))
    fields = [sample_id, len(t_dict), expr_qs, length_qs, num_exon_qs]
    print >>stats_fh, '\t'.join(map(str, fields))


def aggregate(samples, ref_gtf_file, gtf_expr_attr, tmp_dir,
              output_gtf_file, stats_file):
    '''
    Aggregate/merge individual sample GTF files
    '''
    # setup output files
    tmp_file = os.path.join(tmp_dir, 'transcripts.unsorted.gtf')
    tmp_fileh = open(tmp_file, 'w')
    stats_fileh = open(stats_file, 'w')
    # stats file has header
    fields = ['sample_id', 'num_transfrags', 'expr_quantiles',
              'length_quantiles', 'num_exon_quantiles']
    print >>stats_fileh, '\t'.join(fields)
    # aggregate ref gtf
    if ref_gtf_file is not None:
        sample = Sample(ref_gtf_file, Sample.REF_ID)
        sample._id = Sample.REF_ID
        logging.debug('Reference: %s' % ref_gtf_file)
        _aggregate_gtf(sample.gtf_file, sample._id, gtf_expr_attr,
                       tmp_fileh, stats_fileh, is_ref=True)

    # aggregate sample gtfs
    for sample in samples:
        logging.debug('Sample: %s %s' % (sample._id, sample.gtf_file))
        _aggregate_gtf(sample.gtf_file, sample._id, gtf_expr_attr,
                       tmp_fileh, stats_fileh)
    tmp_fileh.close()
    stats_fileh.close()

    # sort merged gtf
    logging.info("Sorting GTF")
    retcode = sort_gtf(tmp_file, output_gtf_file, tmp_dir=tmp_dir)
    if retcode != 0:
        logging.error("Error sorting GTF")
        if os.path.exists(output_gtf_file):
            os.remove(output_gtf_file)
        raise TacoError('Error sorting GTF')
    os.remove(tmp_file)
