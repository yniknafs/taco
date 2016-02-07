import argparse
import logging
import sys

import pysam


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('--frac', type=float, default=0.0)
    parser.add_argument('gtf_file')
    args = parser.parse_args()

    all_t_ids = set()
    t_ids = set()
    for f in pysam.tabix_iterator(open(args.gtf_file), pysam.asGTF()):
        if f.feature == 'transcript':
            t_id = f.transcript_id
            frac = float(f.frac)
            keep = (frac >= args.frac)
            all_t_ids.add(t_id)
            if keep:
                t_ids.add(t_id)
                print str(f)
        elif f.feature == 'exon':
            t_id = f.transcript_id
            assert t_id in all_t_ids
            if t_id in t_ids:
                print str(f)

if __name__ == '__main__':
    sys.exit(main())
