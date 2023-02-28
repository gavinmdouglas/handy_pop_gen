#!/usr/bin/python3

from functions.iupac import check_ambig_match_dict
import sys

codon_to_aa = {"TTT":"F", "TTC":"F",
                "TTA":"L", "TTG":"L",
                "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                "TAT":"Y", "TAC":"Y",
                "TAA":"STOP", "TAG":"STOP",
                "TGT":"C", "TGC":"C",
                "TGA":"STOP",
                "TGG":"W",
                "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                "CAT":"H", "CAC":"H",
                "CAA":"Q", "CAG":"Q",
                "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                "ATT":"I", "ATC":"I", "ATA":"I",
                "ATG":"M",
                "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                "AAT":"N", "AAC":"N",
                "AAA":"K", "AAG":"K",
                "AGT":"S", "AGC":"S",
                "AGA":"R", "AGG":"R",
                "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                "GAT":"D", "GAC":"D",
                "GAA":"E", "GAG":"E",
                "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


def start_codon_present_M_only(seq):
    '''Check that start codon is present in an input DNA sequence of gene at the very beginning.
    If not present will check to see when the first start codon occurs. If any
    ambiguous bases present then will only call it not present if none of
    the possible codons match. Also, return the position of the first start codon if
    present. If there is no canonical start codon, but a non-canonical start codon is present,
    then return the percent of the overall gene that is truncated.'''

    start_codon_position = None

    expected_start_codon = seq[:3]

    if not check_ambig_match_dict(expected_start_codon, codon_to_aa, 'M'):

        exp_start_codon_present = False

        for i in range(3, len(seq), 3):

            codon = seq[i:i + 3]

            if len(codon) != 3:
                continue

            if check_ambig_match_dict(codon, codon_to_aa, 'M'):
                start_codon_position = i

    else:
        exp_start_codon_present = True
        start_codon_position = 0

    if not exp_start_codon_present and start_codon_position:
        leading_percent_truncated = (start_codon_position / len(seq)) * 100
    else:
        leading_percent_truncated = float('NaN')

    return((exp_start_codon_present, start_codon_position, leading_percent_truncated))


def stop_codon_premature_present(seq, start_codon_position):
    '''Check whether premature stop codon is present in an input DNA sequence of
    gene. Also check if the expected stop codon at the end is missing.
    If there are ambiguous bases in a codon
    then premature stop codons must be matched by all possible ambiguous codons,
    while the expected codons need only be matched by one. Also, takes in position
    of start codon, so that if the start codon is missing and/or a few codons in
    that a stop codon isn't called before it.'''

    if len(seq) % 3 != 0:
        print('Warning: gene sequence is not perfectly divisible by three.',
              file = sys.stderr)

    if not start_codon_position:
        start_codon_position = 0

    premature_stop_codon_position = float('NaN')
    percent_truncated = float('NaN')

    stop_codon_calls = set()

    codon_positions = []

    for i in range(start_codon_position, len(seq), 3):

        codon = seq[i:i + 3]

        if len(codon) != 3:
            continue

        codon_positions.append(i)

        if check_ambig_match_dict(codon, codon_to_aa, 'STOP',
                                  all_possible_match = True):
            stop_codon_calls.add(i)

    
    # Re-check tail of gene while allowing for more ambiguous matches to stop
    # codon
    expected_stop_pos = codon_positions[-1]

    final_stop_present = False

    obs_final_codon = seq[expected_stop_pos:expected_stop_pos + 3]

    if len(obs_final_codon) == 3 and check_ambig_match_dict(obs_final_codon, codon_to_aa, 'STOP', all_possible_match = False):

        final_stop_present = True
        stop_codon_calls.remove(expected_stop_pos)


    if len(stop_codon_calls) > 0:
        premature_stop_codon_position = min(list(stop_codon_calls))

        # Twos in numerator just added for clarity - the twos in just general are included as the positions are from the *start* of the codon.
        expected_length = expected_stop_pos + 2 - start_codon_position
        observed_length = premature_stop_codon_position - 2 - start_codon_position

        percent_truncated = ((expected_length - observed_length) / expected_length) * 100

    return((final_stop_present, premature_stop_codon_position, expected_stop_pos, percent_truncated))

