#!/usr/bin/env python3

"""
Script used to remove non-biallelic snps from SYNC files
"""

import os
import pandas as pd
import pathlib
from pathlib import Path
import datetime
import argparse
import numpy as np
from random import sample
import gzip
import sys

# Function to check if a given site is biallelic (i.e., contains only two alleles).
# If the site has more than two alleles in any population, it flags it as non-biallelic.
def check_if_biallelic(sync_line):
    total_count = 0
    count_extra_allele = 0
    count_extra_allele_multiple_pops = 0
    count_extra_allele_perpop = np.array(list([0]*34))  # Array for allele counts across pops (34 pops assumed)
    max_num_popsWith_extra = 1
    max_num_counts_per_pop = 1
    max_proportion_reads_extra_allele = 0.00  # Threshold for extra allele proportion
    line_split = sync_line.strip().split("\t")  # Split the line into columns
    just_pops = line_split[3:]  # Extract population allele data (columns after the first 3)
    allele_counts = [":".join(x.split(':')[:4]) for x in just_pops]  # Only keep the first 4 allele counts
    
    # Check if there's a third allele in any population (biallelic check)
    check_for_third_allele = [np.sort(list(map(int, allele_counts[i].split(':'))))[1] != 0 for i in range(0, len(allele_counts))]
    
    if True in check_for_third_allele:  # If any population has a third allele (i.e., not biallelic)
        check_list = []  # List to store the proportion of extra alleles
        for i in range(0, len(allele_counts)):
            sorted_allele_counts = np.sort(list(map(int, allele_counts[i].split(':'))))  # Sort allele counts
            if sorted_allele_counts[1] != 0:
                add = sorted_allele_counts[1] / (sorted_allele_counts[2] + sorted_allele_counts[3])  # Calculate proportion of extra allele
                check_list.append(round(add, 3))  # Store the proportion rounded to 3 decimal places
            else:
                check_list.append(0)  # If there's no extra allele, append 0
        
        # Count how many populations have extra alleles
        num_pops_NUMALT_3or4 = len([1 for x in check_list if x > 0])
        which_pops = [index + 1 for (index, item) in enumerate(check_list) if item != 0]  # List of pops with extra alleles

        count_extra_allele_perpop = count_extra_allele_perpop + check_list  # Update the allele count per population
        total_count += 1
        count_extra_allele += 1
        
        if num_pops_NUMALT_3or4 > 1:
            count_extra_allele_multiple_pops += 1  # Increment if multiple populations have extra alleles
        
        return True  # Return True to indicate the site is not biallelic
    else:
        total_count += 1
        return False  # Return False if the site is biallelic (i.e., no third allele)

# Main function to process the sync file and filter based on biallelic test.
def main(sync_file, sync_outfile):
    with gzip.open(sync_file, 'rt') as data, gzip.open(sync_outfile, 'wb') as file_to_write:
        for sync_line in data:
            line_split = sync_line.strip().split("\t")  # Split the line into columns
            just_pops = line_split[3:]  # Extract population data
            allele_counts = [":".join(x.split(':')[:6]) for x in just_pops]  # Take first 6 allele counts from each pop
            
            # Check if the line is biallelic
            multiple_allele_test = check_if_biallelic(sync_line)
            
            if multiple_allele_test == False:  # If site is biallelic, write to output
                file_to_write.write(sync_line.encode())  # Write the line to the output file (gzip encoded)
            else:
                pass  # If not biallelic, skip writing the line (do nothing)

# Main execution block, triggered when the script is run directly
if __name__ == "__main__":
    # Get the input and output file paths from command line arguments
    sync_file = sys.argv[1]
    sync_outfile = sync_file[:-4] + "ALT.sync.gz"  # Output file name with "_ALT" suffix
    
    # Call the main function to process the sync file
    main(sync_file, sync_outfile)
