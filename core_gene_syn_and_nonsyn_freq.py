#!/usr/bin/python3

import argparse
import sys
import os
from functions.io_utils import read_fasta
from functions.codon import codon_to_aa
import glob
import json
import itertools
import random
import math
from collections import defaultdict
import pprint


def main():

    parser = argparse.ArgumentParser(

        description='Read in panaroo gene breakdown table and folder containing codon-aligned FASTAs of (single-copy) core genes per species. '
                    'Identify all non-synonymous and synonymous sites and compute the folded site-frequency spectra for each site type separately '
                    '(note that all core genes will be concatenated per strain/genome). '
                    'Codons with > 1 polymorphic sites (across any genome) will be ignored and individual sites will be ignored if they have more than two bases. '
                    'Set number of genomes to subsample to for the specified number of replicates (mean frequencies across replicates will be taken) .'
                    'Note that fewer than the specified subsampling replicates is possible, depending on the possible combination. '
                    'You can avoid subsampling entirely if you set --num_genomes to be the actual number of genomes. '
                    'Also, note that invariant codons will be ignored.',

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
                        help='Path to folder containing files breaking down ',
                        required=False,
                        default=None)

    parser.add_argument('--panaroo',
                        metavar='TABLE',
                        type=str,
                        help='Panaroo output table (gene_presence_absence.csv).',
                        required=True)

    parser.add_argument('--num_genomes',
                        metavar='COUNT',
                        type=int,
                        help='Number of genomes to subsample to before computing site-frequency spectrum (based on set number of replicates)',
                        required=True)

    parser.add_argument('--num_reps',
                        metavar='COUNT',
                        type=int,
                        help='Maximum number of subsampling replicates.',
                        required=True)

    parser.add_argument('--fasta_suffix',
                        metavar='TABLE',
                        type=str,
                        help='File suffix to add to gene id to get FASTA filenames.',
                        required=False,
                        default='.msa.fna')

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

    # The codon info will be re-coded as 0's and 1's per sample
    # (artbitrarily to distinguish the two alleles), to quickly
    # be able to parse the frequencies, which will be separated
    # by syn and nonsyn sites.
    binary_syn_codons = []
    binary_nonsyn_codons = []

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
                clean_seqs[genome_match] = [raw_seq[i: i + 3] for i in range(0, seq_length, 3)]

                if genome_match not in represented_genomes:
                    represented_genomes.add(genome_match)
                else:
                    sys.exit('Error: multiple gene ids matched to same genome for this gene family: ' + gf)

        # Sanity checks on full set of genes and genomes for this gene family:
        if genome_ids != sorted(list(represented_genomes)):
            #print('Skipping this gene family as not all genomes were represented in the FASTA MSA: ' + gf,
            #      file = sys.stderr)
            continue

        if exp_genes_present != len(families_to_genes[gf]):
            sys.exit('Error: not all expected genes found for this gene family: ' + gf)

        # Now loop through every codon position and re-code as appropriate if there is varaition across the samples.
        # (or ignore otherwise if there is no variation, or if there are more than two codons observed).
        num_raw_codons = clean_seqs[genome_ids[0]]

        for codon_i in range(len(num_raw_codons)):
            variant_highest_code = -1
            variant_map = {}
            codon_code = []
            for genome_id in genome_ids:
                codon_instance = clean_seqs[genome_id][codon_i]
                if codon_instance not in variant_map.keys():
                    variant_highest_code += 1
                    variant_map[codon_instance] = variant_highest_code
                codon_code.append(variant_map[codon_instance])

            if len(variant_map) != 2:
                continue
            else:
                codons = list(variant_map.keys())

                skip_codon = False
                for c in codons:
                    if c not in codon_to_aa.keys():
                        skip_codon = True
                if skip_codon:
                    continue

                if codon_to_aa[codons[0]] == codon_to_aa[codons[1]]:
                    binary_syn_codons.append(codon_code)
                else:
                    binary_nonsyn_codons.append(codon_code)

    # After reading in all nonsyn and syn site codes, then 
    # generate all possible combinations of pairs
    all_genome_combos = list(itertools.combinations(list(range(len(genome_ids))), args.num_genomes))

    # Pick up to 100 random pairs (without replacement) from the list of combinations
    random_combos = random.sample(all_genome_combos, min(args.num_reps, len(all_genome_combos)))

    syn_site_freqs = defaultdict(int)
    nonsyn_site_freqs = defaultdict(int)

    for combo_i in random_combos:
        for syn_codon_i in range(len(binary_syn_codons)):
            syn_site_freqs[sum([binary_syn_codons[syn_codon_i][i] for i in combo_i])] += 1

        for nonsyn_codon_i in range(len(binary_nonsyn_codons)):
            nonsyn_site_freqs[sum([binary_nonsyn_codons[nonsyn_codon_i][i] for i in combo_i])] += 1

    max_folded_freq = math.floor(args.num_genomes / 2)

    print('freq\tsyn\tnonsyn')
    for i in range(1, max_folded_freq + 1, 1):
        if i != (args.num_genomes - i):
            syn_mean = (syn_site_freqs[i] + syn_site_freqs[args.num_genomes - i]) / args.num_reps
            nonsyn_mean = (nonsyn_site_freqs[i] + nonsyn_site_freqs[args.num_genomes - i]) / args.num_reps
        else:
            syn_mean = (syn_site_freqs[i]) / args.num_reps
            nonsyn_mean = (nonsyn_site_freqs[i]) / args.num_reps
        print('\t'.join([str(i), str(syn_mean), str(nonsyn_mean)]))


if __name__ == '__main__':
    main()
