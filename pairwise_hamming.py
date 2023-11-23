#!/usr/bin/python3

import argparse
import sys
from functions.io_utils import read_fasta
import itertools
import collections

def main():

    parser = argparse.ArgumentParser(

        description='Compute Hamming distance between each pair of sequences in a FASTA. '
                    'Ignore sites that are gaps (-) in both sequences, but treats a gap only in one sequence as a difference. '
                    'Ignores sites where a base besides A, C, G, T, (or a gap) is present. '
                    'Percent identity returned is just ',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='Path to FASTA file.', required=True)

    # Optionally remove header to make it easy to concatenate output files.
    parser.add_argument('-n', '--no_header', action='store_true',
                        help='Do not write a header in output')

    parser.add_argument('-e', '--extra_annot', type=str,
                        required=False, default=None,
                        help='Value of optional first column to include (e.g., gene name).')

    parser.add_argument('-m', '--id_map', type=str,
                        required=False, default=None,
                        help='Optional tab-delimited file of (1) all FASTA header '
                             'IDs and (2) the corresponding IDs to use in output. '
                             'One ID pair per line.')

    args = parser.parse_args()

    # Read in input sequences.
    seqs = read_fasta(args.input, convert_upper=True, cut_header=True)

    # Convert FASTA IDs if mapfile specified.
    if args.id_map is not None:
        id_map = {}
        with open(args.id_map, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                id_map[line[0]] = line[1]

        previous_id_counts = collections.defaultdict(int)
        tmp = {}

        for seq_id in seqs.keys():
            if seq_id not in id_map.keys():
                sys.exit('Error: ID ' + seq_id + ' not found in mapfile.')

            if previous_id_counts[id_map[seq_id]] > 0:
                new_id = id_map[seq_id] + '_' + str(previous_id_counts[id_map[seq_id]])
                tmp[new_id] = seqs[seq_id]
            else:
                tmp[id_map[seq_id]] = seqs[seq_id]

            previous_id_counts[id_map[seq_id]] += 1

        seqs = tmp

    all_ids = sorted(list(seqs.keys()))

    # Check all sequences are the same length.
    seq_length = len(seqs[all_ids[0]])
    for seq_id in all_ids[1:]:
        if len(seqs[seq_id]) != seq_length:
            sys.exit('Error: All sequences must be the same length.')

    pairwise_combos = list(itertools.combinations(all_ids, 2))

    bases = ['A', 'C', 'G', 'T']

    headerline = 'seq1\tseq2\ttotal_length\tcomparable_length\thamming_dist\tpercent_identity'

    if args.extra_annot is not None:
        headerline =  'annot\t' + headerline
        initial_outline = [args.extra_annot]
    else:
        initial_outline = []

    if not args.no_header:
        print(headerline)

    for combo in pairwise_combos:
        seq1 = seqs[combo[0]]
        seq2 = seqs[combo[1]]

        outline = initial_outline.copy()
        outline += [combo[0], combo[1], str(seq_length)]

        comparable_length = 0
        hamming_dist = 0

        for i in range(seq_length):
            if seq1[i] == '-' and seq2[i] == '-':
                continue
            elif seq1[i] == '-' or seq2[i] == '-':
                comparable_length += 1
                hamming_dist += 1
            elif seq1[i] not in bases or seq2[i] not in bases:
                continue
            elif seq1[i] == seq2[i]:
                comparable_length += 1
            else:
                comparable_length += 1
                hamming_dist += 1

        percent_identity = (1 - (hamming_dist / comparable_length)) * 100

        outline += [str(comparable_length),
                    str(hamming_dist),
                    str(percent_identity)]

        print('\t'.join(outline))


if __name__ == '__main__':
    main()
