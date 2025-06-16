#!/usr/bin/env python3

"""
This script takes a VCF file of SNPs and another VCF file of INDELs. It then
filters the SNP VCF for any SNPs flanking indels. It is different than typical INDEL
filters by taking into account the length of the INDEL as estimated by Freebayes. 
"""


import sys
import argparse

def process_indel_vcf(indel_vcf_content, flanking_bases=10):
    """
    Processes a VCF file containing indels and generates a list of regions 
    that are within a given number of bases of an indel.

    Parameters:
    - indel_vcf_content: A list of strings, each representing a line from the VCF.
    - flanking_bases: The number of bases to include before and after the indel.
    
    Returns:
    - A set of 'chrom_pos' strings representing the regions to mask.
    """
    regions = set()  # Using a set for fast lookup

    for line in indel_vcf_content:
        # Skip header lines
        if line.startswith("#"):
            continue
        
        columns = line.strip().split("\t")
        
        indel_pos = int(columns[1])  # Indel position
        allele_1_len = len(columns[3])  # Length of the reference allele
        
        # Get lengths of alternative alleles (split if multiple alternatives)
        alt_alleles = columns[4].split(",")
        alt_allele_lengths = [len(alt) for alt in alt_alleles]
        max_alt_len = max(alt_allele_lengths, default=0)
        
        # Determine the region to mask
        max_len = max(allele_1_len, max_alt_len)
        start = indel_pos - flanking_bases
        stop = indel_pos + max_len + flanking_bases
        
        # Add positions to the regions set
        for pos in range(start, stop + 1):
            regions.add(f"{columns[0]}_{pos}")

    return regions

def filter_snp_vcf(vcf_file, regions_set, output_file):
    """
    Filters the SNP VCF based on the indel regions and writes the result to a new VCF file.
    
    Parameters:
    - vcf_file: Path to the SNP VCF file.
    - regions_set: A set containing 'chrom_pos' strings representing the regions to mask.
    - output_file: Path to the output filtered VCF file.
    """
    with open(vcf_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write VCF headers
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                break
        
        # Process the SNP VCF content
        for line in infile:
            columns = line.strip().split("\t")
            chrom_pos = f"{columns[0]}_{columns[1]}"  # Combine chrom and position
            
            # Only write to the output if the position is NOT in the regions set
            if chrom_pos not in regions_set:
                outfile.write(line)

    print(f"Filtered VCF saved to: {output_file}")

def main():
    # Read indel VCF content from stdin or a file
    indel_vcf_file = sys.argv[1]  # Indel VCF file
    snp_vcf_file = sys.argv[2]    # SNP VCF file
    output_vcf_file = sys.argv[3] # Output VCF file
    output_vcf_file_suffix = (f"{output_vcf_file}_IND")

    with open(indel_vcf_file, 'r') as infile:
        indel_vcf_content = infile.readlines()

    # Process indel VCF and generate regions
    regions_set = process_indel_vcf(indel_vcf_content)

    # Filter SNP VCF based on indel regions and write the result
    filter_snp_vcf(snp_vcf_file, regions_set, output_vcf_file_suffix)

if __name__ == "__main__":
    main()


#To run this 
""" 
python script.py \
sard_ncbi_all_Qual21_indels.vcf \
sard_ncbi_all_Qual30_biSnps_STB_HAPs.vcf \
sard_ncbi_all_Qual30_biSnps_STB_HAPs_IND.vcf
"""