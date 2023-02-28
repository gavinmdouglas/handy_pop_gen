#!/usr/bin/python3

import argparse
import sys
import os

from functions.io_utils import read_fasta


def main():

    parser = argparse.ArgumentParser(

    description="Reads in FASTA and write out two specified sequences to a new file in \"axt\" format (which is used by KaKs_Calculator).",

    formatter_class=argparse.RawDescriptionHelpFormatter)


    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        required=True, help="Path to FASTA file")

    parser.add_argument("-n", "--name", metavar="NAME", type=str,
                        required=True, help="Name to print at top of axt outfile.")

    parser.add_argument("--ids2keep", metavar="FILE", type=str, required=True,
                        help="File with ids to keep - one per line.")

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        required=True, help="Path to output axt file")

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header=True)
    
    ids2keep = []
    with open(args.ids2keep, 'r') as ids_file:
        for line in ids_file:
            
            line_split = line.split()

            if len(line_split) > 0:
                ids2keep.append(line_split[0])

    if len(ids2keep) != 2:
        sys.exit('Stopping - exactly two sequence ids to retain need to be specified.')

    for id2keep in ids2keep:
        if id2keep not in seqs.keys():
            sys.exit('Stopping, id ' + id2keep + ' not found in input FASTA.')

    out_axt = open(args.output, "w")
    out_axt.write(args.name + "\n")

    for id2keep in ids2keep:
        out_axt.write(seqs[id2keep] + "\n")

    out_axt.close()


if __name__ == '__main__':
    main()

