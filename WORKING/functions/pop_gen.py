#!/usr/bin/python3

from math import sqrt, comb
import numpy as np
import itertools
import sys
from collections import Counter
import codon

def tajimas_d_and_diversity_metrics(in_seqs):
    '''Return n, S, Watterson's theta, nucleotide diversity, and Tajima's D as a tuple
    for a list of input sequences.'''

    (theta_pi, S) = calc_theta_pi_and_num_segregating(in_seqs)

    N = len(in_seqs)

    (tajimas_d, wattersons_theta) = tajimas_d_and_wattersons_theta(theta_pi, S, N)

    return((N, S, wattersons_theta, theta_pi, tajimas_d))


def diversity_syn_vs_nonsyn(in_seqs):
    '''Calculate and return Theta Pi as well as the number of segregating
    sites as a tuple, for non-syn, any syn, and four-fold syn sites separately.
    Note: this function will not count polymorphisms where pairwise seqs do
    not have bases that are on of A, C, G, or T.'''

    bases = ['A', 'C', 'G', 'T']

    # Dereplicate sequences
    inseqs_breakdown = Counter(in_seqs)

    if len(inseqs_breakdown.keys()) == 1:
        return((0, 0))

    pairwise_comparisons = list(itertools.combinations(inseqs_breakdown.keys(), 2))

    nonsyn_segregating_sites = set()
    nonsyn_pairwise_diffs = []

    anysyn_segregating_sites = set()
    anysyn_pairwise_diffs = []

    fourfoldsyn_segregating_sites = set()
    fourfoldsyn_pairwise_diffs = []

    for combo in pairwise_comparisons:

        nonsyn_num_diff = 0
        anysyn_num_diff = 0
        fourfoldsyn_num_diff = 0

        seq1 = combo[0]
        seq2 = combo[1]

        if len(seq1) != len(seq2):
            print(seq1)
            print(seq2)
            sys.exit('Stopping - lengths of input sequences do not match. They should be aligned!')

        if len(seq1) % 3 != 0:
            sys.exit('Stopping - alignment length not a multiple of three: should be a codon alignment!')

        for i in range(3, len(seq), 3):

            codon1 = seq1[i:(i + 3)]
            codon2 = seq2[i:(i + 3)]

            

            if seq1[i] != seq2[i]:
                if seq1[i] in bases and seq2[i] in bases:
                    num_diff += 1
                    segregating_sites.add(i)

        # Add the num diff for every time the non-dereplicated sequences would have been compared.
        for compare_i in range(inseqs_breakdown[seq1] * inseqs_breakdown[seq2]):
            pairwise_diffs.append(num_diff)

    # Also add in num diff's of 0 for redundant sequences.
    for unique_seq in inseqs_breakdown.keys():
        unique_seq_freq = inseqs_breakdown[unique_seq]
        
        if unique_seq_freq == 1:
            continue

        for duplicated_seq_combo in range(comb(unique_seq_freq, 2)):
            nonsyn_pairwise_diffs.append(0)
            anysyn_pairwise_diffs.append(0)
            fourfoldsyn_pairwise_diffs.append(0)

    return((np.mean(pairwise_diffs), len(segregating_sites)))


def calc_theta_pi_and_num_segregating(in_seqs):
    '''Calculate and return Theta Pi as well as the number of segregating
    sites as a tuple. Note: this function will not count polymorphisms
    where pairwise seqs do not have bases that are on of A, C, G, or T.'''

    bases = ['A', 'C', 'G', 'T']

    # Dereplicate sequences
    inseqs_breakdown = Counter(in_seqs)

    if len(inseqs_breakdown.keys()) == 1:
        return((0, 0))

    pairwise_comparisons = list(itertools.combinations(inseqs_breakdown.keys(), 2))

    segregating_sites = set()

    pairwise_diffs = []

    for combo in pairwise_comparisons:

        num_diff = 0

        seq1 = combo[0]
        seq2 = combo[1]

        if len(seq1) != len(seq2):
            print(seq1)
            print(seq2)
            sys.exit('Stopping - lengths of input sequences do not match. They should be aligned!')

        for i in range(len(seq1)):
        
            if seq1[i] != seq2[i]:
                if seq1[i] in bases and seq2[i] in bases:
                    num_diff += 1
                    segregating_sites.add(i)

        # Add the num diff for every time the non-dereplicated sequences would have been compared.
        for compare_i in range(inseqs_breakdown[seq1] * inseqs_breakdown[seq2]):
            pairwise_diffs.append(num_diff)

    # Also add in num diff's of 0 for redundant sequences.
    for unique_seq in inseqs_breakdown.keys():
        unique_seq_freq = inseqs_breakdown[unique_seq]
        
        if unique_seq_freq == 1:
            continue

        for duplicated_seq_combo in range(comb(unique_seq_freq, 2)):
            pairwise_diffs.append(0)

    return((np.mean(pairwise_diffs), len(segregating_sites)))


def num_pairwise_diff_bases(bases):
    '''Will return the number of differences for a set of bases. Note: this
    function assumes that all bases are all one of A, C, T, G.'''

    num_diff = 0

    for possible_base in ['A', 'C', 'T', 'G']:
        
        num_specific_base = sum(x == possible_base for x in bases)
        num_other_base = len(bases) - num_specific_base

        if num_specific_base == 0 or num_other_base == 0:
            continue

        num_diff += num_specific_base * num_other_base
        bases = [b for b in bases if b != possible_base]

    return(num_diff)


def tajimas_d_and_wattersons_theta(theta_pi, num_seg_sites, n):
    '''Compute Tajima's D based on theta pi, number of segregating sites, and
    num_samples (n). Will return Watterson's theta as the second returned
    element as well.'''

    a1 = sum(1 / x for x in range(1, n))

    a2 = sum(1 / (x**2) for x in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))

    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))

    c1 = b1 - (1 / a1)

    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1**2)  

    e1 = c1 / a1

    e2 = c2 / (a1**2 + a2)

    expected_sd = sqrt(e1 * num_seg_sites + e2 * num_seg_sites * (num_seg_sites - 1))

    wattersons_theta = num_seg_sites / a1

    if expected_sd == 0:
        tajimas_d = float("NaN")
    else:
        tajimas_d = (theta_pi - wattersons_theta) / expected_sd

    return(tajimas_d, wattersons_theta)

