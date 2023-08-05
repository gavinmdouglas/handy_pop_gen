#!/usr/bin/python3

import argparse
from collections import defaultdict
import math
import sys

from functions.codon import codon_to_aa
from functions.io_utils import float_prop_arg, read_fasta
from functions.nucleotide import check_DNA_RNA_IUPAC, tally_gc_content


def main():

    parser = argparse.ArgumentParser(

        description='''
        Computes the Codon Adaptation Index of Species (CAIS), based on the description in Weibel et al.
        2023 (https://doi.org/10.7554/eLife.87335.1).This index encompasses both per-species GC content
        and differences in amino acid frequencies.

        Requires an input FASTA file of all protein-coding sequences (separately) for a given species.
        The species' amino acid frequencies will be parsed from this input (along with the codon
        frequencies per amino acid). Note that start codons will be trimmed from genes unless
        --assume_M_start is specified, which differs from the implementation in Weibel et al. 2023. This
        is to avoid the issue that multiple start codons are used by prokaryotes, which were the intended
        use cases for this script. Note that the GC content should be specified with the --gc argument
        in most cases (e.g., to input the genome-wide GC content), but if not specified, the GC content
        inferred from the input genes will be used.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to FASTA file')

    parser.add_argument('-a', '--aa_freq_table', metavar='TABLE', type=str,
                        required=True,
                        help='Path to tab-delimited table with '
                             'amino acid frequencies (as proportions) based '
                             'on all species being compared.')

    parser.add_argument('--gc', metavar='FLOAT', type=float_prop_arg,
                        required=False, default=None,
                        help='Optional GC-content (as proportion) for species. '
                             'If not provided, will be calculated from input FASTA, '
                             'which will not be representative of the genome-wide GC content '
                             'as used in Weibel et al. 2023.')

    parser.add_argument('--assume_M_start', action='store_true',
                        help='Assume that all sequences should start with methionine codon. '
                             'Otherwise start codons will be ignored, and will not be '
                             'considered when computing CAIS.')

    parser.add_argument('--print_gene_warnings', action='store_true',
                        help='Print warnings about genes that are not divisible by three in '
                             'length and that (optionally) do not begin with a methionine codon.')

    parser.add_argument('-n', '--name', metavar='NAME', type=str,
                        required=False, default=None,
                        help='Optional name to output with CAIS (to make it easy to concatenate '
                             'with other outputs). If unspecified, the input FASTA file name '
                             'will be output.')

    args = parser.parse_args()

    if not args.name:
        import os

    seqs = read_fasta(args.fasta, cut_header=True)

    if len(seqs) == 0:
        sys.exit('Stopping - at least one gene sequence must be input.')

    if len(seqs) < 500:
        print('Read in ' + str(len(seqs)) + ' sequences. '
              'This is a relatively small number of sequences: '
              'make sure that you input a FASTA containing the '
              'protein-coding DNA sequences for this species.',
              file=sys.stderr)

    # Loop through all input sequences and perform the following steps:
    # 1. Check that length is divisible by three, otherwise ignore that gene and throw warning.
    # 2. Convert to uppercase.
    # 3. Check that all bases are IUPAC DNA/RNA characters.
    # 4. Tally up GC content (unless GC content argument specified).
    # 5. Tally up frequency of all codons (and for each amino acid they encode).
    # Note that the start codons will be ignored, unless the --assume_M_start argument was provided.

    gc_tally = 0
    total_length = 0
    codon_freq = defaultdict(int)
    aa_freq = defaultdict(int)
    total_codons = 0
    num_ignored_seqs = 0

    for seq_id, seq in seqs.items():

        if len(seqs[seq_id]) % 3 != 0:
            if args.print_gene_warnings:
                print('Warning: Sequence ' + seq_id + ' is not divisible by three and will be ignored.',
                      file=sys.stderr)
            num_ignored_seqs += 1
            continue

        seqs[seq_id] = seq.upper()

        # Replace any U characters with T.
        seqs[seq_id] = seq.replace('U', 'T')

        if args.assume_M_start:
            if seq[0:3] != 'ATG':
                if args.print_gene_warnings:
                    print('Warning: Sequence ' + seq_id + ' does not start with ATG '
                          '(which is expected as --assume_M_start was specified) and will be ignored.',
                          file=sys.stderr)
                num_ignored_seqs += 1
                continue
        else:
            # Trim start codon.
            seqs[seq_id] = seqs[seq_id][3:]

        check_DNA_RNA_IUPAC(seqs[seq_id])

        if not args.gc:
            gc_tally += tally_gc_content(seqs[seq_id], check_IUPAC=False, type='count')
            total_length += len(seqs[seq_id])

        # Loop over all codons in sequence and increment tallies.
        for i in range(0, len(seqs[seq_id]), 3):
            codon = seqs[seq_id][i:i + 3]

            # Skip if contains any ambiguous characters.
            if codon not in codon_to_aa.keys():
                continue

            codon_freq[codon] += 1
            aa_freq[codon_to_aa[codon]] += 1
            total_codons += 1

    if args.assume_M_start:
        print(str(num_ignored_seqs) + ' sequences were ignored due to being '
              'non-divisible by three and/or not beginning with ATG.',
              file=sys.stderr)

    elif num_ignored_seqs > 0:
        print(str(num_ignored_seqs) + ' sequences were ignored due to being '
              'non-divisible by three.', file=sys.stderr)

    if not args.gc:
        gc_content = gc_tally / total_length
    else:
        gc_content = args.gc

    if gc_content == 0 or gc_content == 1:
        sys.exit('Stopping - GC content is ' + str(gc_content) + '. '
                 'Non-zero GC and AT content is required to calculate CAIS.')

    at_content = 1 - gc_content

    # Check for any missing codons and give them a frequency of one count.
    for codon in codon_to_aa.keys():
        if codon_freq[codon] == 0:
            codon_freq[codon] = 1
            total_codons += 1
            aa_freq[codon_to_aa[codon]] += 1

    # Convert aa and codon counts to proportions.
    for aa in aa_freq.keys():
        aa_freq[aa] = aa_freq[aa] / total_codons

    for codon in codon_freq.keys():
        codon_freq[codon] = codon_freq[codon] / total_codons

    # Read in amino acid frequencies from table.
    all_aa = set('ACDEFGHIKLMNPQRSTVWY*')
    aa_cross_species_freq = dict()
    total_aa_freq = 0
    with open(args.aa_freq_table, 'r') as aa_filehandle:
        for line in aa_filehandle:
            line = line.rstrip()
            line = line.split('\t')
            aa_code = line[0]
            in_freq = float(line[1])
            total_aa_freq += in_freq

            if aa_code not in all_aa:
                sys.exit('Stopping - amino acid code ' + aa_code + ' is not valid. '
                         'Valid amino acid codes are: ' + ', '.join(all_aa))

            if aa_code in aa_cross_species_freq.keys():
                sys.exit('Stopping - amino acid code ' + aa_code +
                         ' is duplicated in table.')

            aa_cross_species_freq[line[0]] = float(line[1])

    if abs(1 - total_aa_freq) > 0.00001:
        sys.exit('Stopping - total amino acid frequency is ' + str(total_aa_freq) +
                 ', but should be 1.')

    # Compute relative synonymous codon usage (RSCU) for each codon.
    # First, get the probability of observing each codon based on the GC content.
    raw_codon_prob_by_gc = dict()
    for codon in codon_to_aa.keys():
        codon_gc_count = tally_gc_content(codon, check_IUPAC=False, type='count')
        raw_codon_prob_by_gc[codon] = (gc_content ** codon_gc_count) * (at_content ** (3 - codon_gc_count))

    # Then get the expected probability of observing each codon relative to all other codons for an AA.
    # To do so, first need to get mapping for AA to codons.
    aa_to_codon = defaultdict(set)
    for codon, aa in codon_to_aa.items():
        aa_to_codon[aa].add(codon)

    # Then compute the relative synonymous codon usage.
    codon_to_RSCU = dict()
    for aa in aa_to_codon:
        aa_summed_raw_prob = 0
        aa_summed_obs_freq = 0

        for codon in aa_to_codon[aa]:
            aa_summed_raw_prob += raw_codon_prob_by_gc[codon]
            aa_summed_obs_freq += codon_freq[codon]

        for codon in aa_to_codon[aa]:

            codon_exp_by_gc = raw_codon_prob_by_gc[codon] / aa_summed_raw_prob
            codon_obs_by_gc = codon_freq[codon] / aa_summed_obs_freq
            codon_to_RSCU[codon] = codon_obs_by_gc / codon_exp_by_gc

    # Next, compute the species-specific amino acid weighting factor.
    aa_weights = dict()
    for aa in aa_cross_species_freq.keys():
        aa_weights[aa] = aa_cross_species_freq[aa] / aa_freq[aa]

    # Last, compute the CAIS score.
    CAIS = 0
    for codon in codon_to_RSCU:
        CAIS += codon_freq[codon] * aa_weights[codon_to_aa[codon]] * math.log(codon_to_RSCU[codon])

    CAIS = math.exp(CAIS)

    if args.name:
        name = args.name
    else:
        name = os.path.basename(name)

    print(name + '\t' + str(CAIS))


if __name__ == '__main__':
    main()
