#!/usr/bin/python3

import egglib
import sys

# Get pi estimates at synonymous sites only and other relevant stats. 
# Only argument should be input FASTA.
def main():

    aln = egglib.io.from_fasta(sys.argv[1], alphabet=egglib.alphabets.DNA)

    aln_codons = egglib.tools.to_codons(aln)

    coding_diversity_out = egglib.stats.CodingDiversity(align = aln_codons)

    num_stop = coding_diversity_out.num_codons_stop

    struct = egglib.struct_from_labels(aln_codons)
    cs = egglib.stats.ComputeStats()
    cs.set_structure(struct)
    cs.add_stats('Pi', 'lseff', 'nseff')
    stats = cs.process_sites(coding_diversity_out.sites_S)

    print('\t'.join([sys.argv[1], str(stats['Pi']), str(stats['lseff']), str(stats['nseff'])]))

    stats = cs.process_sites(coding_diversity_out.sites_NS)

    print('\t'.join([sys.argv[1], str(stats['Pi']), str(stats['lseff']), str(stats['nseff'])]))


if __name__ == '__main__':

    main()
