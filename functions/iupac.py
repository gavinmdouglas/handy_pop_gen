#!/usr/bin/python3

import itertools
import sys

iupac_bases = dict()

iupac_bases['M'] = ['A', 'C']
iupac_bases['R'] = ['A', 'G']
iupac_bases['W'] = ['A', 'T']

iupac_bases['S'] = ['C', 'G']
iupac_bases['Y'] = ['C', 'T']
iupac_bases['K'] = ['G', 'T']

iupac_bases['V'] = ['A', 'C', 'G']
iupac_bases['H'] = ['A', 'C', 'T']
iupac_bases['D'] = ['A', 'G', 'T']

iupac_bases['B'] = ['C', 'G', 'T']

iupac_bases['N'] = ['A', 'C', 'G', 'T']


def check_ambig_match_dict(input_string, dict2check, expected_val, all_possible_match=False):
    '''Given an input string with DNA characters, will check if it is the input
    dictionary or not. If it is then it will check whether the expected value is
    returned or not. Otherwise, assuming that ambiguous DNA characters are
    present, it will generate all  possible strings made up of the four non-ambiguous bases. 
    Will then check to see if any of these non-ambigous strings match as keys to the input 
    dictionary and if so if the expected value for that key is found. The
    purpose of this function that was originally in mind was to check to see
    if codons that match the expected value were detected. By default will
    return True if *any* of the possible non-ambiguous strings match the
    expected value. When all_possible_match=True it will only return True when
    all possible non-ambiguous strings match the expected value.'''

    # First check if input string is in dictionary. Return whether it matches
    # expected value or not if so.
    if input_string in dict2check:
        return(dict2check[input_string] == expected_val)

    possible_bases = []
    ambig_base_count = 0

    for i in range(len(input_string)):
        base_i = input_string[i]
        if base_i in iupac_bases:
            possible_bases.append(iupac_bases[base_i])
            ambig_base_count += 1
        else:
            possible_bases.append([base_i])

    if ambig_base_count == 0:
        print(input_string, file = sys.stderr)
        sys.exit('Error - input string is not found in input dictionary and '
                 'also does not contain any ambiguous bases.')

    non_ambig_combos = list(itertools.product(*possible_bases))

    match_count = 0
    for non_ambig in non_ambig_combos:
        non_ambig = ''.join(non_ambig)

        if dict2check[non_ambig] == expected_val:
            match_count += 1

    if not all_possible_match and match_count > 0:
        return(True)
    elif all_possible_match and match_count == len(non_ambig_combos):
        return(True)
    else:
        return(False)

