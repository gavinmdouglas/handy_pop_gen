#!/usr/bin/python3

import argparse
from functions.io_utils import read_fasta
from functions.codon import codon_to_aa
from collections import defaultdict
import sys


def main():

    parser = argparse.ArgumentParser(

        description='''
        Read in FASTA of coding sequences (one gene per entry) and
        return proportions of all AAs based on translated codons
        (including stop codons).
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-f', '--fasta',
                        metavar='FILE',
                        type=str,
                        help='Path to input FASTA file with DNA sequence(s).',
                        required=True)

    parser.add_argument('--ignore_start_codon', action='store_true',
                        help='Ignore start codon when computing frequencies '
                             '(it will be trimmed from each sequence).')

    args = parser.parse_args()

    in_seqs = read_fasta(args.fasta,
                         convert_upper=True)

    total_codons = 0
    aa_counts = defaultdict(int)

    for seq_id, nucl in in_seqs.items():

        # Remove any gap characters.
        nucl = nucl.replace('.', '').replace('-', '')

        if len(nucl) % 3 != 0:
            print('Skipping sequence ' + seq_id +
                  'as it is not divisible by three.',
                  file=sys.stderr)
            continue

        if args.ignore_start_codon:
            # Trim start codon.
            nucl = nucl[3:]

        for i in range(0, len(nucl), 3):
            codon = nucl[i:i + 3]

            # Skip if contains any ambiguous characters.
            if codon not in codon_to_aa.keys():
                continue
            else:
                aa_counts[codon_to_aa[codon]] += 1
                total_codons += 1

    print('aa\tprop')
    for aa in sorted(aa_counts.keys()):
        print(aa + '\t' + str(aa_counts[aa] / total_codons))


if __name__ == '__main__':
    main()
