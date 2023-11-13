#!/usr/bin/python3

import argparse
import sys
from functions.io_utils import read_fasta, write_fasta
import os

def main():

    parser = argparse.ArgumentParser(

        description='Read in folder containing subfolders of codon-aligned FASTAs of (single-copy) core genes per species. '
                    'Each subfolder is assumed to be the actual gene family ID. '
                    'Create combined multiple-sequence alignment for all single-copy, core genes. '
                    'Each genome is expected to be present in each gene family alignment once (i.e., they are assumed to be dereplicated already!).',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-s', '--species_folder',
                        metavar='FOLDER',
                        type=str,
                        help='Path to folder containing subfolders named for each gene family.',
                        required=True)

    parser.add_argument('--fasta_suffix',
                        metavar='STRING',
                        type=str,
                        help='File suffix to add to gene ID to get FASTA filenames.',
                        required=False,
                        default='.msa.fna')

    parser.add_argument('--genome_id_fields',
                        metavar='INT',
                        type=int,
                        help='Number of fields in FASTA file headers after splitting by _ to determine genome ID.',
                        required=False,
                        default=3)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output FASTA.',
                        required=True)

    args = parser.parse_args()

   # Get lists of all subfolders in directory.
    subfolders = [f.path for f in os.scandir(args.species_folder) if f.is_dir()]

    core_genome_alignment = dict()
    genome_ids = []
    first_subfolder = True
    for subfolder in subfolders:
        gf = os.path.basename(subfolder)

        fasta_file = subfolder + '/' + gf + args.fasta_suffix

        if not os.path.isfile(fasta_file):
            sys.exit('Error: FASTA file not found: ' + fasta_file)

        raw_seqs = read_fasta(fasta_file)

        skip = False
        uniq_seq_ids = list(raw_seqs.keys())
        seq_length = len(raw_seqs[uniq_seq_ids[0]])
        for uniq_seq_id in uniq_seq_ids[1:]:
            if seq_length != len(raw_seqs[uniq_seq_id]):
                print('Error: sequences for this gene are of different length, but they should be aligned: ' + gf, file = sys.stderr)
                skip = True

        if skip:
            continue

        if seq_length % 3 != 0:
            print('Error: codon-aligned sequences for this gene are not divisible by 3: ' + gf, file = sys.stderr)
            continue

        if len(raw_seqs) == 1:
            continue

        # Loop through sequences, remove "_N" from ids and add in any identical seqs
        # (use genome id as key, as raw gene id is not relevant).
        # Also run sanity check that all expected genes and genomes are represented.
        represented_genomes = set()

        # Sometimes weirdly a gene specified as identical in the mapfile will be separate
        # in the actual final msa, so need to make sure the same gene id was not already
        # looped through (this happened at least once).
        past_derep = set()
        clean_seqs = dict()
        for raw_seq_id, raw_seq in raw_seqs.items():
            genome_id = '_'.join(raw_seq_id.split('_')[:args.genome_id_fields])

            if genome_id in past_derep:
                sys.exit('Error: duplicate genome ID found in this gene family: ' + gf)
            else:
                past_derep.add(genome_id)

            clean_seqs[genome_id] = raw_seq

            if genome_id not in represented_genomes:
                represented_genomes.add(genome_id)
            else:
                sys.exit('Error: multiple gene ids matched to same genome for this gene family: ' + gf)
    
        if first_subfolder:
            genome_ids = sorted(list(clean_seqs.keys()))       
            first_subfolder = False

            for genome_id in genome_ids:
                core_genome_alignment[genome_id] = clean_seqs[genome_id]
        else:
            # Sanity checks on full set of genes and genomes for this gene family:
            if genome_ids != sorted(list(represented_genomes)):
                # Skip as for some reason the genome ids do not match!
                print('Warning: mismatch in observed and expected genomes for this gene family: ' + gf, file = sys.stderr)
        
            # Then add each sequence to full core genome alignment.
            for genome_id in genome_ids:
                core_genome_alignment[genome_id] += clean_seqs[genome_id]

    # Final sanity check that all sequences are the same length.
    seq_length = len(core_genome_alignment[genome_ids[0]])
    for genome_id in genome_ids[1:]:
        if seq_length != len(core_genome_alignment[genome_id]):
            sys.exit('Error: not all final core genome alignments are the same length!')

    if seq_length % 3 != 0:
        sys.exit('Error: final core genome alignment not divisible by three.')

    write_fasta(core_genome_alignment, args.output)


if __name__ == '__main__':
    main()
