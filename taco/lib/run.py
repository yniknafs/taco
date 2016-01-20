'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
@author: yniknafs
'''
import os
import argparse
import logging
import pickle
import json

from base import Sample
from aggregate import aggregate
from assemble import assemble

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2015"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class Args:
    VERBOSE = False
    GUIDED = False
    GTF_EXPR_ATTR = 'TPM'
    MIN_FRAG_LENGTH = 200
    MAX_FRAG_LENGTH = 400
    FRAC_ISOFORM = 0.01
    MAX_ISOFORMS = 100
    OUTPUT_DIR = 'assemblyline'
    RESUME = False
    PROG = 'taco'
    DESCRIPTION = 'Meta-assembly of RNA-Seq datasets'

    @staticmethod
    def create():
        parser = argparse.ArgumentParser(description=Args.DESCRIPTION)
        parser.add_argument('-v', '--verbose', dest='verbose',
                            action="store_true",
                            default=Args.VERBOSE,
                            help='Enabled detailed logging '
                            '(for debugging)')
        parser.add_argument('--gtf-expr-attr',
                            dest='gtf_expr_attr',
                            default=Args.GTF_EXPR_ATTR,
                            metavar='ATTR',
                            help='GTF attribute field containing '
                            'expression [default=%(default)s]')
        parser.add_argument('--min-frag-length', dest='min_frag_length',
                            type=int, metavar='N',
                            default=Args.MIN_FRAG_LENGTH,
                            help='Length (bp) of smallest valid fragment '
                            'across all experiments [default=%(default)s]')
        parser.add_argument('--max-frag-length', dest='max_frag_length',
                            type=int, metavar='N',
                            default=Args.MAX_FRAG_LENGTH,
                            help='Length (bp) of largest valid fragment '
                            'across all experiments [default=%(default)s]')
        parser.add_argument('--ref-gtf', dest='ref_gtf_file',
                            metavar='<gtf_file>',
                            default=None,
                            help='Option reference GTF file of "true" '
                            'validated transcripts to facilitate guided '
                            'assembly and/or noise filtering')
        parser.add_argument('--guided', dest='guided', action='store_true',
                            default=Args.GUIDED,
                            help='Enable guided assembly (requires a '
                            'reference GTF to be specified using '
                            '--ref-gtf)')
        parser.add_argument('--frac-isoform',
                            dest='frac_isoform', type=float,
                            metavar='X',
                            default=Args.FRAC_ISOFORM,
                            help='Report transcript isoforms with '
                            'expression fraction >=X (0.0-1.0) of the total '
                            'gene expression [default=%(default)s]')
        parser.add_argument('--max-isoforms',
                            dest='max_isoforms', type=int, metavar='N',
                            default=Args.MAX_ISOFORMS,
                            help='Maximum isoforms to report for each '
                            'gene [default=%(default)s]')
        parser.add_argument("-o", "--output-dir", dest="output_dir",
                            metavar='DIR',
                            default=Args.OUTPUT_DIR,
                            help='directory where output files will be '
                            'stored (if already exists then --resume must '
                            'be specified) [default=%(default)s]')
        parser.add_argument('--resume', dest='resume',
                            action='store_true',
                            default=Args.RESUME,
                            help='resume a previous run in <dir>')
        parser.add_argument('sample_file', nargs='?')
        return parser

    @staticmethod
    def load(filename):
        return pickle.load(open(filename))

    @staticmethod
    def dump(args, filename):
        pickle.dump(args, open(filename, 'wb'))

    @staticmethod
    def log(args, func=logging.info):
        logging.info('%s version %s' % (Args.PROG, __version__))
        spacer = '-' * 78
        fmt = '{:<35}{:<35}'
        func(spacer)
        func(fmt.format('verbose logging:', str(args.verbose)))
        func(fmt.format('output directory:', str(args.output_dir)))
        func(fmt.format('reference GTF file:', str(args.ref_gtf_file)))
        func(fmt.format('guided mode:', str(args.guided)))
        func(fmt.format('GTF expression attribute:', args.gtf_expr_attr))
        func(fmt.format('min fragment length:', args.min_frag_length))
        func(fmt.format('max fragment length:', args.max_frag_length))
        func(fmt.format('fraction isoform:', args.frac_isoform))
        func(fmt.format('max isoforms:', args.max_isoforms))
        return

    @staticmethod
    def parse():
        parser = Args.create()
        args = parser.parse_args()
        # check if we are trying to resume a previous run
        if args.resume:
            if not os.path.exists(args.output_dir):
                parser.error("Output directory '%s' does not exist" %
                             args.output_dir)
            # check that basic files exist in directory
            results = Results(args.output_dir)
            can_resume = os.path.exists(results.args_file)
            can_resume = can_resume and os.path.exists(results.status_file)
            can_resume = can_resume and os.path.exists(results.sample_file)
            if not can_resume:
                parser.error("Output directory '%s' not valid" %
                             args.output_dir)
        else:
            if os.path.exists(args.output_dir):
                parser.error("Output directory '%s' already exists" %
                             args.output_dir)
            if args.sample_file is None:
                parser.error("sample file not specified")
            if not os.path.exists(args.sample_file):
                parser.error("sample file %s not found" % (args.sample_file))
            args.sample_file = os.path.abspath(args.sample_file)

            if args.min_frag_length <= 0:
                parser.error("min_transcript_length <= 0")
            if args.max_frag_length <= 0:
                parser.error("max_frag_length <= 0")
            if (args.frac_isoform <= 0) or (args.frac_isoform > 1):
                parser.error("frac_isoform out of range (0.0-1.0)")
            if (args.max_isoforms < 1):
                parser.error("max_isoforms <= 0")

            if args.ref_gtf_file is not None:
                if not os.path.exists(args.ref_gtf_file):
                    parser.error("reference GTF file %s not found" %
                                 (args.ref_gtf_file))
                args.ref_gtf_file = os.path.abspath(args.ref_gtf_file)
            elif args.guided:
                parser.error('Guided assembly mode (--guided) requires '
                             'reference GTF (--ref-gtf)')
        return args


class Results(object):
    TMP_DIR = 'tmp'
    STATUS_FILE = 'status.json'
    ARGS_FILE = 'args.pickle'
    SAMPLE_FILE = 'samples.txt'
    TRANSFRAGS_GTF_FILE = 'transfrags.gtf'
    TRANSFRAGS_FAIL_GTF_FILE = 'transfrags.fail.gtf'
    AGGREGATE_STATS_FILE = 'aggregate_stats.txt'

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.tmp_dir = \
            os.path.join(output_dir, Results.TMP_DIR)
        self.args_file = \
            os.path.join(output_dir, Results.ARGS_FILE)
        self.status_file = \
            os.path.join(output_dir, Results.STATUS_FILE)
        self.sample_file = \
            os.path.join(output_dir, Results.SAMPLE_FILE)
        self.transfrags_gtf_file = \
            os.path.join(output_dir,
                         Results.TRANSFRAGS_GTF_FILE)
        self.transfrags_fail_gtf_file = \
            os.path.join(output_dir,
                         Results.TRANSFRAGS_FAIL_GTF_FILE)
        self.aggregate_stats_file = \
            os.path.join(output_dir,
                         Results.AGGREGATE_STATS_FILE)


class Status(object):
    FIELDS = ('create', 'aggregate', 'assemble')

    def __init__(self):
        for f in Status.FIELDS:
            setattr(self, f, False)

    def write(self, filename):
        d = dict((f, getattr(self, f)) for f in Status.FIELDS)
        json.dump(d, open(filename, 'w'))

    @staticmethod
    def load(filename):
        d = json.load(open(filename))
        self = Status()
        for f in Status.FIELDS:
            setattr(self, f, d[f])
        return self


class Run(object):

    @staticmethod
    def create():
        self = Run()
        # parse command line args
        args = Args.parse()
        if args.resume:
            self.results = Results(args.output_dir)
            self.status = Status.load(self.results.status_file)
            self.samples = Sample.parse_tsv(self.results.sample_file,
                                            header=True)
            self.args = Args.load(self.results.args_file)
        else:
            self.args = args
            self.results = Results(args.output_dir)
            self.status = Status()
            self.samples = Sample.parse_tsv(self.args.sample_file)

        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s - %(levelname)s - %(message)s")

        # create output directories
        results = self.results
        if not os.path.exists(results.output_dir):
            logging.debug("Creating output directory '%s'" %
                          (results.output_dir))
            os.makedirs(results.output_dir)
        if not os.path.exists(results.tmp_dir):
            logging.debug("Creating tmp directory '%s'" % (results.tmp_dir))
            os.makedirs(results.tmp_dir)

        if not args.resume:
            # write command line args
            Args.dump(self.args, self.results.args_file)
            # write samples
            Sample.write_tsv(self.samples, self.results.sample_file)
            # update status and write to file
            self.status.create = True
            self.status.write(self.results.status_file)
        return self

    def log_args(self):
        Args.log(self.args)

    def aggregate(self):
        '''
        Aggregate/merge individual sample GTF files
        '''
        r = self.results
        a = self.args
        samples = self.samples

        aggregate(samples,
                  ref_gtf_file=a.ref_gtf_file,
                  gtf_expr_attr=a.gtf_expr_attr,
                  tmp_dir=r.tmp_dir,
                  output_gtf_file=r.transfrags_gtf_file,
                  stats_file=r.aggregate_stats_file)

        # update status and write to file
        self.status.aggregate = True
        self.status.write(self.results.status_file)

    def assemble(self):
        assemble(gtf_file=self.results.transfrags_gtf_file,
                 output_dir=self.args.output_dir)
        # update status and write to file
        # self.status.assemble = True
        # self.status.write(self.results.status_file)
