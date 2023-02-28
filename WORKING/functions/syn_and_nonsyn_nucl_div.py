#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.dn_ds import pairwise_dnds
import itertools
import pandas as pd
import numpy as np

def main():

    parser = argparse.ArgumentParser(

        description='Read in FASTA file of codon-aligned CDS sequences and compute pi, theta, and Tajima\'s D (the last three are for all sites, syn, four-fold syn, and non-syn sites, respectively).',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='Path to codon-aligned FASTA.', required=True)

    parser.add_argument('--include_header', action='store_true',
                        help='Include a header for the output.',
                        required=False)

    args = parser.parse_args()


    # Read in input sequences.
    seqs = read_fasta(args.input)

    if len(seqs.keys()) == 0:
        sys.exit("Stopping: no sequences were found in the input file.")
    elif len(seqs.keys()) == 1:
        sys.exit("Stopping: there is only one input sequence present.")

    
    
    if args.include_header:
        print('\t'.join(['dN', 'dS', 'dN/dS']))




if __name__ == '__main__':
    main()
