#!/usr/bin/python3

import argparse
from functions.io_utils import read_fasta
from functions.nucleotide import tally_gc_content


def main():

    parser = argparse.ArgumentParser(

        description='''
        Read in FASTA and return GC content as
        proportion based on all input DNA sequences.
        Note that IUPAC ambiguous characters are allowed
        and will be counted as fractional GC (e.g., N
        will be counted as a 0.5 GC count rather than 0).
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-f', '--fasta',
                        metavar='FILE',
                        type=str,
                        help='Path to input FASTA file with DNA sequence(s).',
                        required=True)

    args = parser.parse_args()

    in_seqs = read_fasta(args.fasta,
                         convert_upper=True)

    total_length = 0
    total_gc = 0

    for nucl in in_seqs.values():

        # Remove any gap characters.
        nucl = nucl.replace('.', '').replace('-', '')

        total_length += len(nucl)

        total_gc += tally_gc_content(nucl,
                                     check_IUPAC=True,
                                     type='count')

    print(total_gc / total_length)


if __name__ == '__main__':
    main()
