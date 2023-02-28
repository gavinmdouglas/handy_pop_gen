#!/usr/bin/python3

import argparse
import sys
import os
from functions.io_utils import read_fasta
from functions.dn_ds import pairwise_dnds
import itertools
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

        description='Read in FASTA file of CDS sequences and compute dN, dS, and dN/dS for '
                    'each pairwise comparison. Will return mean of these pairwise comparisons '
                    '(and optionally a table of all individual pairwise comparisons as well). '
                    'NOTE: This approach does not correct for mutational biases (e.g., different '
                    'rates of transition vs transversion mutations), as implemented in Li 1993.',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='Path to codon-aligned FASTA.', required=True)

    parser.add_argument('--min_num_subs_total', metavar='COUNT', type=int,
                        help='Min number of subs (combined syn. and non-syn.) needed before including the pairwise comparison.',
                        required=False, default=0)

    parser.add_argument('--min_num_subs_each', metavar='COUNT', type=int,
                        help='Min number of subs (for each of syn. and non-syn. subs, i.e., they both need to have x subs) needed before including the pairwise comparison.',
                        required=False, default=0)

    parser.add_argument('--output_pairwise', metavar='OUTPUT', type=str,
                        help='Optional name of output file to write values for each pairwise comparison. '
                             'These will not be output unless this option is specified.',
                        required=False)

    parser.add_argument('--exclude_header', action='store_true',
                        help='Exclude header when writing out the mean values.',
                        required=False)

    args = parser.parse_args()

    # Read in input sequences.
    seqs = read_fasta(args.input, convert_upper=True)

    if len(seqs.keys()) == 0:
        sys.exit("Stopping: no sequences were found in the input file.")
    elif len(seqs.keys()) == 1:
        sys.exit("Stopping: there is only one input sequence present.")

    # Get mean dn, ds, and dn/ds based on all pairwise comparisons.
    pairwise_combos = list(itertools.combinations(seqs, 2))

    dnds = pd.DataFrame(columns=['n_subs', 'n_sites',
                                 's_subs', 's_sites',
                                 'dn', 'ds', 'dnds'],
                        index=pd.MultiIndex.from_tuples(pairwise_combos,
                                                        names=('seq1', 'seq2')))
    for combo in pairwise_combos:
        dnds.loc[combo, :] = list(pairwise_dnds(seqs[combo[0]],
                                                seqs[combo[1]]))

    dnds = dnds.loc[(dnds.n_subs + dnds.s_subs >= args.min_num_subs_total), :]
    dnds = dnds.loc[(dnds.n_subs >= args.min_num_subs_each), :]
    dnds = dnds.loc[(dnds.s_subs >= args.min_num_subs_each), :]

    mean_n_subs = np.nanmean(np.array(dnds.loc[:, 'n_subs']), dtype='float32')
    mean_s_subs = np.nanmean(np.array(dnds.loc[:, 's_subs']), dtype='float32')
    mean_n_sites = np.nanmean(np.array(dnds.loc[:, 'n_sites']), dtype='float32')
    mean_s_sites = np.nanmean(np.array(dnds.loc[:, 's_sites']), dtype='float32')
    mean_dn = np.nanmean(np.array(dnds.loc[:, 'dn']), dtype='float32')
    mean_ds = np.nanmean(np.array(dnds.loc[:, 'ds']), dtype='float32')
    mean_dnds = np.nanmean(np.array(dnds.loc[:, 'dnds']), dtype='float32')

    if not args.exclude_header:
        print('\t'.join(['infile', 'mean_n_subs', 'mean_n_sites',
                         'mean_s_subs', 'mean_s_sites',
                         'mean_dn', 'mean_ds', 'mean_dnds']))

    mean_output = [os.path.basename(args.input), mean_n_subs, mean_n_sites,
                   mean_s_subs, mean_s_sites,
                   mean_dn, mean_ds, mean_dnds]
    mean_output = [str(val) for val in mean_output]
    print('\t'.join(mean_output))

    if args.output_pairwise:
        dnds.to_csv(args.output_pairwise,
                    sep='\t', header=True, na_rep='NA')


if __name__ == '__main__':
    main()
