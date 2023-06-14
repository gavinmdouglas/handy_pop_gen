#!/usr/bin/python3

import argparse
import sys

from functions.io_utils import read_fasta


def main():

    parser = argparse.ArgumentParser(

        description="Read in FASTA of aligned sequences. "
                    "Print out table with rows as each position in the "
                    "alignment and columns as each separate sequence id. "
                    "Values will be the actual bases in each sequence at the "
                    "specified position. This script is intended to be used to "
                    "run sanity checks on alignments to clearly visualize "
                    "differences at specific positions. "
                    "Can specify whether all positions or only variable positions "
                    "be written out. First position output is 1-based (i.e., does not start at 0).",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to FASTA file')

    parser.add_argument('--variable_only', action='store_true',
                        help='Only output sites with at least two different bases.',
                        required=False)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header=True)

    if len(seqs) <= 1:
        sys.exit('Stopping - at least two aligned sequences must be output.')

    seq_ids = sorted(seqs.keys())
    obs_length = len(seqs[seq_ids[0]])
    for seq_id in seq_ids[1:]:
        if obs_length != len(seqs[seq_id]):
            sys.exit('Error - sequence lengths differ, but aligned sequences expected.')

    print(','.join(['position'] + seq_ids))

    pos = 0
    for i in range(obs_length):
        pos += 1
        diff_bases = set()
        outline = []
        for seq_id in seq_ids:
            outline.append(seqs[seq_id][i])
            diff_bases.add(seqs[seq_id][i])

        if args.variable_only and len(set(outline)) == 1:
            continue

        print(','.join([str(pos)] + outline))


if __name__ == '__main__':
    main()
