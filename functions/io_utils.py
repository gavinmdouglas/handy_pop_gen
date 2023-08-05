#!/usr/bin/python3

import argparse
import gzip
import sys
import textwrap

def read_ids(filename):
    '''Read ids from a file into a list (one id per line)'''
    ids = list()
    with open(filename, 'r') as id_file:
        for line in id_file:
            ids.append(line.rstrip())

    return(ids)

def read_fasta(filename, cut_header=False, convert_upper=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    def parse_fasta_lines(file_handle):
        for line in file_handle:

            line = line.rstrip()

            if len(line) == 0:
                continue

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Make sure that sequence id is not already in dictionary.
                if name in seq:
                    sys.stderr("Stopping due to duplicated id in file: " + name)

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        with gzip.open(filename, "rt") as fasta_in:
            parse_fasta_lines(fasta_in)
    else:
        with open(filename, "r") as fasta_in:
            parse_fasta_lines(fasta_in)

    if convert_upper:
        for seq_name in seq.keys():
            seq[seq_name] = seq[seq_name].upper()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")

    out_fasta.close()


def read_fastq_headers(filepath):
        
        # Capture every 4th line of FASTQ (after first line)
        lineno = 4
        header_lines = []

        with open(filepath, "r") as fastq_in:
            
            for line in fastq_in:
                
                if lineno == 4:
                    header_lines.append(line.rstrip())
                    lineno = 1
                
                else:
                    lineno += 1

        return(header_lines)


def read_vcf_variant_bases(in_vcf, only_poly=False):
    '''Read VCF and return dictionary with coordinate as key and all 
    observed bases at polymorphic / non-reference site as values. Set
    only_poly=True to only return sites where there are two alleles
    segregating (i.e., fixed non-ref sites will not be returned).'''

    variant_sites = {}

    with open(in_vcf, 'r') as vcffile:
        for vcf_line in vcffile:

            if vcf_line[0] == "#":
                continue

            vcf_line = vcf_line.rstrip()
            vcf_split = vcf_line.split()

            contig_name = vcf_split[0]
            position = int(vcf_split[1]) - 1

            ref_base = vcf_split[3]
            alt_bases = vcf_split[4].split(",")

            genotypes = vcf_split[9].split("/")

            observed_bases = []

            if '0' in genotypes:
                observed_bases.append(ref_base)

            for alt_i in range(1, 4):
                if str(alt_i) in genotypes:
                    observed_bases.append(alt_bases[alt_i - 1])

            if only_poly and len(observed_bases) == 1:
                continue

            variant_sites[contig_name + "|" + str(position)] = observed_bases

    return(variant_sites)

def float_prop_arg(arg):
    try:
        value = float(arg)
        if 0.0 <= value <= 1.0:
            return value
        raise argparse.ArgumentTypeError(f"{arg} is not in the range [0.0, 1.0]")
    except ValueError:
        raise argparse.ArgumentTypeError(f"{arg} is not a valid float value")
