#!/usr/bin/python3

import argparse
import sys
import os
import pysam
import math
from collections import defaultdict

def n_choose_k(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def site_within_sample_pi(allele_depth):
    '''Takes in depth of all alleles present at a site,
    and computes the proportion of pairwise read comparisons that differ.
    '''
    # Return 0 if there is only one allele at this site.
    if len(allele_depth) == 1:
        return 0
    elif len(allele_depth) == 0:
        return None

    # Compute total number of pairwise comparisons.
    # And number of comparisons where the reads differ.
    total_comparisons = 0
    total_diff = 0


    for i in range(len(allele_depth)):
        # Keep track of number of pairwise comparisons that are between reads with the same allele.
        total_comparisons += n_choose_k(allele_depth[i], 2)

        for j in range(i + 1, len(allele_depth)):
            num_diff = allele_depth[i] * allele_depth[j]
            total_comparisons += num_diff
            total_diff += num_diff
    
    # Compute proportion of pairwise comparisons that differ.
    return(total_diff / total_comparisons)


def site_between_sample_pi(alleles1, allele_depth1,
                           alleles2, allele_depth2):
    '''Similar function as site_within_sample_pi above, but takes
    in alleles and matching allele depth for two samples and computes
    the mean pairwise divergence.'''

    # Return None if either sample has 0 alleles at this site.
    if len(allele_depth1) == 0 or len(allele_depth2) == 0:
        return None
    
    # Make sure length of each list matches the number of alleles.
    assert len(alleles1) == len(allele_depth1)
    assert len(alleles2) == len(allele_depth2)

    # Keep track of total number of comparisons and number of comparisons that differ.
    total_comparisons = 0
    total_diff = 0

    # Loop through all alleles for sample 1.
    for i in range(len(alleles1)):
        allele1 = alleles1[i]
        allele1_depth = allele_depth1[i]

        # Loop through all alleles for sample 2.
        for j in range(len(alleles2)):
            allele2 = alleles2[j]
            allele2_depth = allele_depth2[j]

            # Compute number of pairwise comparisons and number of comparisons that differ.
            num_comparisons = allele1_depth * allele2_depth
            total_comparisons += num_comparisons

            if allele1 != allele2:
                total_diff += num_comparisons

    return(total_diff / total_comparisons)

def main():

    # Read in command line arguments.
    parser = argparse.ArgumentParser(
        description='Compute within and between sample nucleotide diversity for a VCF file.')
    
    parser.add_argument('--bcf', type=str, required=True,
                        help='Path to BCF file.')
    
    parser.add_argument('--min_total_depth', type=int, default=20,
                        help='Minimum total depth to consider a site.')
    
    parser.add_argument('--min_alt_depth', type=int, default=5,
                        help='Minimum alternate depth to consider a site.')
    
    args = parser.parse_args()

    # Read in BCF file with pysam.
    bcf_in = pysam.VariantFile(args.bcf)

    # Get mapping of raw to clean sample IDs.
    sample_map = dict()
    for sample in bcf_in.header.samples:
        clean_sample = os.path.basename(sample).split('.')[0]
        sample_map[sample] = clean_sample

    # Get list of contigs in BCF file.
    all_contigs = bcf_in.header.contigs

    print('contig', 'comparison_type', 'sample1', 'sample2', 'pi', 'num_seg', 'num_sites', sep='\t')

    # Loop through contigs and then loop through all sites per contig.
    for contig_id in all_contigs:
        
        within_sample_site_pi = defaultdict(list)
        within_sample_site_nseg = defaultdict(int)

        between_sample_site_pi = defaultdict(list)
        between_sample_site_nseg = defaultdict(int)

        for record in bcf_in.fetch(contig_id):

            # Parse all samples for this site and determine which have sufficient coverage.
            # Note that reads corresponding to an allele with < min_alt_depth are not counted towards total depth.

            sample_info = dict()

            for sample in record.samples:

                sample_GT = record.samples[sample]['GT']
                sample_AD = record.samples[sample]['AD']
                sample_depth = 0
                filt_allele_depth = []
                filt_allele = []

                # Figure out the total number of unique genotypes, sorted numerically.
                sample_GT = sorted(list(set(sample_GT)))
                if None in sample_GT:
                    sample_GT.remove(None)
                if len(sample_GT) == 0:
                    continue
                
                # Figure out depth per alleles and re-compute total depth.
                for i in range(len(sample_AD)):
                    if sample_AD[i] is not None and sample_AD[i] >= args.min_alt_depth:
                        sample_depth += sample_AD[i]
                        filt_allele_depth.append(sample_AD[i])
                        filt_allele.append(i)

                if sample_depth < args.min_total_depth:
                    continue
                else:
                    
                    # Get clean ID.
                    clean_sample = sample_map[sample]

                    sample_info[clean_sample] = {'filt_DP': sample_depth,
                                                'AD': filt_allele_depth,
                                                'alleles': filt_allele}
                    
                    within_sample_site_pi[clean_sample].append(site_within_sample_pi(sample_info[clean_sample]['AD']))

                    if len(sample_info[clean_sample]['alleles']) > 1:
                        within_sample_site_nseg[clean_sample] += 1

            filt_samples = sorted(list(sample_info.keys()))

            if len(filt_samples) > 1:
                for i in range(len(filt_samples)):
                    for j in range(i + 1, len(filt_samples)):
                        between_sample_site_pi[(filt_samples[i], filt_samples[j])].append(site_between_sample_pi(sample_info[filt_samples[i]]['alleles'],
                                                                                                                sample_info[filt_samples[i]]['AD'],
                                                                                                                sample_info[filt_samples[j]]['alleles'],
                                                                                                                sample_info[filt_samples[j]]['AD']))
                        
                        if len(sample_info[filt_samples[i]]['alleles']) > 1 and len(sample_info[filt_samples[j]]['alleles']) > 1 or sample_info[filt_samples[i]]['alleles'][0] != sample_info[filt_samples[j]]['alleles'][0]:
                            between_sample_site_nseg[(filt_samples[i], filt_samples[j])] += 1

        # Then loop through and compute the mean across all sites for each within and between sample comparison.
        for sample in within_sample_site_pi:
            num_sites = len(within_sample_site_pi[sample])
            num_seg = within_sample_site_nseg[sample]
            within_pi = sum(within_sample_site_pi[sample]) / num_sites
            print(contig_id, 'within', sample, sample, str(within_pi), str(num_seg), str(num_sites), sep='\t')

        for sample_pair in between_sample_site_pi:
            num_sites = len(between_sample_site_pi[sample_pair])
            num_seg = between_sample_site_nseg[sample_pair]
            between_pi = sum(between_sample_site_pi[sample_pair]) / num_sites
            print(contig_id, 'between', sample_pair[0], sample_pair[1], str(between_pi), str(num_seg), str(num_sites), sep='\t')

if __name__ == '__main__':
    main()