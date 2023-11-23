#!/usr/bin/python3

import argparse
import sys
from functions.io_utils import read_fasta, write_fasta
import gzip

def main():

    parser = argparse.ArgumentParser(

        description='Read in Panaroo gene breakdown table and folder containing aligned genes per species. '
                    'Create combined multiple-sequence alignment for all single-copy, core genes. '
                    'This differs from other scripts in this repository, as it just assumes all core gene families '
                    'are already aligned and in the same folder, regardless of whether this is a codon alignment or not'
                    '(and with the original gene IDs rather than genome IDs).',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-f', '--fastas',
                        metavar='FOLDER',
                        type=str,
                        help='Path to folder containing codon-aligned FASTAs',
                        required=True)

    parser.add_argument('--panaroo',
                        metavar='TABLE',
                        type=str,
                        help='Gzipped Panaroo output table (gene_presence_absence.csv.gz).',
                        required=True)

    parser.add_argument('-o', '--output_fasta',
                        metavar='OUTPUT',
                        type=str,
                        help='Output FASTA filename.',
                        required=True)

    parser.add_argument('-c', '--core_bed',
                        metavar='OUTPUT',
                        type=str,
                        help='Bedfile to indicate all core gene families, and their coordinates in final alignment.',
                        required=True)

    parser.add_argument('-p', '--prefix',
                        metavar='PREFIX',
                        type=str,
                        help='Prefix to add to gene family IDs to get FASTA IDs.',
                        required=False,
                        default='')

    parser.add_argument('-s', '--suffix',
                        metavar='SUFFIX',
                        type=str,
                        help='Suffix to add to gene family IDs to get FASTA IDs.',
                        required=False,
                        default='.fa')

    args = parser.parse_args()

    # Get mapping of individual gene calls to genome ids.
    # Also keep track of all expected genes per gene family (for single-copy, core gene families)
    genes_to_genome = dict()
    families_to_genes = dict()

    with gzip.open(args.panaroo, 'rt') as panaroo_fh:
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

    # Then read through codon-aligned FASTA.
    # Run sanity checks, by confirming that all sequences are the same size,
    # and that all genomes are represented.

    core_genome_alignment = dict()
    for genome in genome_ids:
        core_genome_alignment[genome] = ''

    gene_start_pos = 0
    bed_outlines = []

    singlecopy_core_genefamilies = sorted(list(families_to_genes.keys()))
    for gf in singlecopy_core_genefamilies:

        fasta_file = args.fastas + '/' + args.prefix + gf + args.suffix

        raw_seqs = read_fasta(fasta_file)

        uniq_seq_ids = list(raw_seqs.keys())

        # Check sequence lengths.
        gene_length = len(raw_seqs[uniq_seq_ids[0]])
        for uniq_seq_id in uniq_seq_ids[1:]:
            if gene_length != len(raw_seqs[uniq_seq_id]):
                sys.exit('Error: sequences for this gene are of different length, but they should be aligned: ' + gf)

        # Loop through and add to full alignment.
        # Also, run sanity check that all expected genes and genomes are represented.
        represented_genomes = set()

        for seq_id, seq in raw_seqs.items():
            genome_id = genes_to_genome[seq_id]

            core_genome_alignment[genome_id] += seq

            represented_genomes.add(genome_id)

        # Sanity checks on full set of genes and genomes for this gene family:
        if genome_ids != sorted(list(represented_genomes)):
            sys.exit('Not all genomes represented for gene ' + gf)

        # Prep for bed output and increment running start position.
        bed_outlines.append('\t'.join(['core', str(gene_start_pos), str(gene_start_pos + gene_length), gf]))
        gene_start_pos += gene_length

    # Final check that all sequences in alignment are of equal length.
    exp_length = len(core_genome_alignment[genome_ids[0]])
    for genome_id in genome_ids[1:]:
        if exp_length != len(core_genome_alignment[genome_id]):
            sys.exit('Error: final aligned sequences are not of equal length!')

    # Write out final FASTA.
    write_fasta(core_genome_alignment, args.output_fasta)

    # Write out BED.
    with open(args.core_bed, 'w') as bed_fh:
        for bed_outline in bed_outlines:
            bed_fh.write(bed_outline + '\n')


if __name__ == '__main__':
    main()
