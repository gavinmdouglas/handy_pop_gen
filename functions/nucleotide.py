#!/usr/bin/python3

import sys


def tally_gc_content(in_seq, check_IUPAC=True, type='proportion'):
    '''
    Compute GC content of sequence
    (as a proportion or count, indicated by the 'type' parameter).
    Will account for IUPAC ambiguity codes to compute partial GC counts.
    '''

    if type != 'count' and type != 'proportion':
        type_error = 'type must be either "count" or "proportion"'
        raise ValueError(type_error)

    # Check that all bases are IUPAC characters.
    if check_IUPAC:
        check_DNA_RNA_IUPAC(in_seq)

    # Initialize GC count.
    tally_gc_content = 0

    # Iterate over bases in sequence.
    for base in in_seq:

        # If base is G, C, or S, increment GC count.
        if base in 'GCSgcs':
            tally_gc_content += 1

        # Skip if non-GC.
        elif base in 'ATUWatuw':
            continue

        # If base is IUPAC ambiguity code that is G/C or
        # a non-GC base, then increment GC count by 0.5.
        elif base in 'RYKMNrykmn':
            tally_gc_content += 0.5

        # Perform similar checks for whether base is 2/3 GC.
        elif base in 'BVbv':
            tally_gc_content += 2/3

        # Or 1/3 GC.
        elif base in 'DHdh':
            tally_gc_content += 1/3

    # Return GC content as count or proportion.
    if type == 'count':
        return tally_gc_content
    elif type == 'proportion':
        return tally_gc_content / len(in_seq)


def check_DNA_RNA_IUPAC(in_seq):
    '''
    Sanity check that all bases match IUPAC characters only.
    Throws error if non-IUPAC characters are found, otherwise returns None.
    '''

    IUPAC_bases = 'ACGTURYSWKMBDHVN'
    non_IUPAC = set()
    for base in in_seq:
        if base not in IUPAC_bases:
            non_IUPAC.add(base)

    if len(non_IUPAC) > 0:
        sys.exit('Stopping - non-IUPAC RNA/DNA characters found in sequence: ' + ' '.join(sorted(list(non_IUPAC))))
    else:
        return None
