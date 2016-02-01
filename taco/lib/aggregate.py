'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer
'''
import os
import logging
import collections
import operator

from base import Sample, TacoError
from gtf import GTF, GTFError, sort_gtf
from stats import scoreatpercentile

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2015"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def _read_gtf(gtf_file, sample_id, gtf_expr_attr, is_ref=False):
    t_id_map = {}
    t_dict = {}
    exon_dict = collections.defaultdict(lambda: [])
    cur_t_id = 1
    for f in GTF.parse(open(gtf_file)):
        if f.feature == 'transcript':
            t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
            if t_id in t_id_map:
                m = 'GTF "%s" transcript_id "%s" not unique' % (gtf_file, t_id)
                raise GTFError(m)
            # rename transcript id
            new_t_id = "%s.T%d" % (sample_id, cur_t_id)
            t_id_map[t_id] = new_t_id
            cur_t_id += 1
            # save attributes
            if is_ref:
                expr = '0.0'
            else:
                expr = f.attrs[gtf_expr_attr]
            attrs = ((GTF.Attr.TRANSCRIPT_ID, new_t_id),
                     (GTF.Attr.SAMPLE_ID, sample_id),
                     (GTF.Attr.REF, str(int(is_ref))),
                     (GTF.Attr.EXPRESSION, expr))
            f.attrs = dict(attrs)
            t_dict[new_t_id] = f
        elif f.feature == 'exon':
            t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
            # lookup new transcript id
            new_t_id = t_id_map[t_id]
            # store exon feature
            f.attrs = {GTF.Attr.TRANSCRIPT_ID: new_t_id}
            exon_dict[new_t_id].append(f)
    return t_dict, exon_dict


def add_sample_gtf(sample, gtf_expr_attr, output_fileh, stats_fileh,
                   is_ref=False):
    '''
    Reads and renames transfrags
    Normalizes expression by total filtered expression
    '''
    # read gtf file into dict of transcripts
    t_dict, exon_dict = \
        _read_gtf(sample.gtf_file, sample._id, gtf_expr_attr, is_ref)

    exprs = []
    lengths = []
    num_exons = []
    for t_id in t_dict.iterkeys():
        t_feature = t_dict[t_id]
        exon_features = exon_dict[t_id]
        # save expression values
        expr = float(t_feature.attrs[GTF.Attr.EXPRESSION])
        exprs.append(expr)
        # save transfrag length and number of exons
        lengths.append(sum((f.end - f.start) for f in exon_features))
        num_exons.append(len(exon_features))
        # write transcript
        print >>output_fileh, str(t_feature)
        # sort features (exons) by start position
        exon_features.sort(key=operator.attrgetter('start'))
        # write exons
        for feature in exon_features:
            print >>output_fileh, str(feature)

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
    fields = [sample._id, len(t_dict), expr_qs, length_qs, num_exon_qs]
    print >>stats_fileh, '\t'.join(map(str, fields))


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
              'length_quantiles']
    print >>stats_fileh, '\t'.join(fields)

    # aggregate ref gtf
    if ref_gtf_file is not None:
        sample = Sample(ref_gtf_file, Sample.REF_ID)
        sample._id = Sample.REF_ID
        logging.debug('Reference: %s' % ref_gtf_file)
        add_sample_gtf(sample, gtf_expr_attr, tmp_fileh, stats_fileh,
                       is_ref=True)

    # aggregate sample gtfs
    for sample in samples:
        logging.debug('Sample: %s %s' % (sample._id, sample.gtf_file))
        add_sample_gtf(sample, gtf_expr_attr, tmp_fileh, stats_fileh)
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
