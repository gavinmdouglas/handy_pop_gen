#!/usr/bin/python3

import argparse
import sys
from functions.io_utils import read_fasta, write_fasta
import json


def main():

    parser = argparse.ArgumentParser(

        description='Read in panaroo gene breakdown table and folder containing HyPhy-prepared codon-aligned FASTAs of (single-copy) core genes per species. '
                    'Create combined multiple-sequence alignment for all single-copy, core genes.',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-f', '--fastas',
                        metavar='FOLDER',
                        type=str,
                        help='Path to folder containing codon-aligned FASTAs',
                        required=True)

    parser.add_argument('--maps',
                        metavar='FOLDER',
                        type=str,
                        help='Path to folder containing files breaking down genes to representative copies.',
                        required=False,
                        default=None)

    parser.add_argument('--panaroo',
                        metavar='TABLE',
                        type=str,
                        help='Panaroo output table (gene_presence_absence.csv).',
                        required=True)

    parser.add_argument('--fasta_suffix',
                        metavar='STRING',
                        type=str,
                        help='File suffix to add to gene id to get FASTA filenames.',
                        required=False,
                        default='.msa.fna')

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Output FASTA filename.',
                        required=True)

    args = parser.parse_args()

    # Get mapping of individual gene calls to genome ids.
    # Also keep track of all expected genes per gene family (for single-copy, core gene families)
    genes_to_genome = dict()
    families_to_genes = dict()

    with open(args.panaroo, 'r') as panaroo_fh:
        panaroo_header = panaroo_fh.readline()
        panaroo_header = panaroo_header.rstrip()
        header_split = panaroo_header.split(',')
        genome_ids = sorted(header_split[3:])
        for panaroo_line in panaroo_fh:

            panaroo_line = panaroo_line.rstrip()

            panaroo_line_split = panaroo_line.split(',')

            if len(panaroo_line_split) != len(header_split):
                sys.exit('Error - line does not contain the expected number of fields:\n' + panaroo_line)

            # Read through once to figure out if this is single-copy,
            # core gene.
            singlecopy_core = True
            for i in range(3, len(panaroo_line_split), 1):
                if panaroo_line_split[i] == '' or ';' in panaroo_line_split[i] or '_len' in panaroo_line_split[i] or '_stop' in panaroo_line_split[i] or '_refound' in panaroo_line_split[i]:
                    singlecopy_core = False
                    break

            if not singlecopy_core:
                continue

            # For single-copy, core genes, parse as described.
            gene_set = set()
            for i in range(3, len(panaroo_line_split), 1):
                gene = panaroo_line_split[i]
                gene_set.add(gene)
                genes_to_genome[gene] = header_split[i]

            gene_family = panaroo_line_split[0]
            families_to_genes[gene_family] = gene_set

    # Loop through each single-copy, core gene family identified
    # and build concatenated core genome alignment per genome.

    # First, read through mapping of gene ids to identical copies
    # for this gene family (in the "maps" folder).
    # (As based on the currently used approach identical sequences
    # are de-replicated in the final MSA, with the number of copies
    # indicated by _N in the gene names)

    # Then read through codon-aligned FASTA.
    # Run sanity checks, by confirming that all sequences are the same size
    # and multiples of 3 in length. Then dereplicate identical sequences,
    # make sure all expected genes are present, and that all expcted
    # genomes encode one of these genes, and then add the codon info.

    core_genome_alignment = dict()
    for genome in genome_ids:
        core_genome_alignment[genome] = ''

    singlecopy_core_genefamilies = sorted(list(families_to_genes.keys()))
    for gf in singlecopy_core_genefamilies:
        json_file = args.maps + '/' + gf + '.fna_copies.json'
        with open(json_file) as f:
            json_raw = f.read()
        gene_copies_raw = json.loads(json_raw)

        # Sometimes for some reason the representative copy id is different in the final MSA compared to the JSON file.
        # For this reason, get map between every member and whole set, just for simplicity.
        gene_copies = {}
        for rep in gene_copies_raw.keys():
            all_members = list(gene_copies_raw[rep].values())
            for member in all_members:
                gene_copies[member] = all_members

        fasta_file = args.fastas + '/' + gf + '.msa.fna'

        raw_seqs = read_fasta(fasta_file)

        uniq_seq_ids = list(raw_seqs.keys())
        seq_length = len(raw_seqs[uniq_seq_ids[0]])
        for uniq_seq_id in uniq_seq_ids[1:]:
            if seq_length != len(raw_seqs[uniq_seq_id]):
                sys.exit('Error: sequences for this gene are of different length, but they should be aligned: ' + gf)

        if seq_length % 3 != 0:
            sys.exit('Error: codon-aligned sequences for this gene are not divisible by 3: ' + gf)

        # Loop through sequences, remove "_N" from ids and add in any identical seqs
        # (use genome id as key, as raw gene id is not relevant).
        # Also run sanity check that all expected genes and genomes are represented.
        clean_seqs = {}
        exp_genes_present = 0
        represented_genomes = set()

        # Sometimes weirdly a gene specified as identical in the mapfile will be separate
        # in the actual final msa, so need to make sure the same gene id was not already
        # looped through (this happened at least once).
        past_derep = set()

        for raw_seq_id, raw_seq in raw_seqs.items():
            clean_id = '_'.join(raw_seq_id.split('_')[:-1])

            for derep_id in gene_copies[clean_id]:

                if derep_id in past_derep:
                    continue
                else:
                    past_derep.add(derep_id)

                if derep_id in families_to_genes[gf]:
                    exp_genes_present += 1
                else:
                    sys.exit('Error: gene (' + derep_id + ') not expected for gene family: ' + gf)

                genome_match = genes_to_genome[derep_id]

                # Transform sequences into list of codons.
                clean_seqs[genome_match] = raw_seq

                if genome_match not in represented_genomes:
                    represented_genomes.add(genome_match)
                else:
                    sys.exit('Error: multiple gene ids matched to same genome for this gene family: ' + gf)

        # Sanity checks on full set of genes and genomes for this gene family:
        if genome_ids != sorted(list(represented_genomes)):
            # Skip as for some reason the genome ids do not match!
            continue

        if exp_genes_present != len(families_to_genes[gf]):
            sys.exit('Error: not all expected genes found for this gene family: ' + gf)

        # Then add each sequence to full core genome alignment.
        for genome_id in genome_ids:
            core_genome_alignment[genome_id] += clean_seqs[genome_id]

    write_fasta(core_genome_alignment, args.output)


if __name__ == '__main__':
    main()
