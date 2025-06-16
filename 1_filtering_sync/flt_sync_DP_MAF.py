#!/usr/bin/env python3

#Script used to filter a SYNC file for read depth and minimum minor allele frequency (MAF)

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

# Function to determine the number of populations and create a default max depth array.
# This is done by reading the first line of the sync file, extracting the number of pops,
# and initializing a default max depth for each pop.
def determine_num_pops(sync_file):
    with gzip.open(sync_file, 'rt') as data:
        first_line = data.readline().strip().split("\t")  # Read first line and split by tabs
        num_pops = len(first_line[3:])  # Count the number of populations (columns beyond the first 3)
        max_depth_default = [100] * num_pops  # Default depth of 100 for each population
        print(len(max_depth_default), max_depth_default)  # Debug: print default depths
    return num_pops, max_depth_default

# Function to extract allele counts and depths from a single line of the sync file.
def get_AC_depths(line, max_depths):
    line_split = line.strip().split("\t")  # Split the line into tabs
    just_pops = line_split[3:]  # Extract population data (from the 4th column onward)
    allele_counts = [":".join(x.split(':')[:4]) for x in just_pops]  # Get first 4 elements from each pop's data
    depths = [np.sum(list(map(int, allele_counts[i].split(':'))))  # Sum the counts for depth calculation
              for i in range(0, len(max_depths))]
    return allele_counts, depths

# Function to test if any population depth is below the minimum required
def test_min_depth(depths, min_depth=10):
    depth_test_min = any([x < min_depth for x in depths])  # Check if any depth is below min_depth
    return depth_test_min

# Function to test if any population depth exceeds the max depth for that population
def test_max_depth(depths, max_depths):
    depth_test_max = (False in [i >= j for i, j in zip(max_depths, depths)])  # Compare each depth with max_depths
    return depth_test_max

# Function to calculate Minor Allele Frequency (MAF) and test against a minimum average threshold
def test_min_maf(allele_counts, max_depths, min_maf_average):
    major_Afs = [sorted(list(map(int, allele_counts[i].split(':'))))[-1]
                 for i in range(0, len(max_depths))]  # Get the major allele frequency
    minor_Afs = [sorted(list(map(int, allele_counts[i].split(':'))))[-2]
                 for i in range(0, len(max_depths))]  # Get the minor allele frequency

    if (0 not in major_Afs):  # Avoid division by zero if major allele is zero
        mafs = [i/j for i in minor_Afs for j in major_Afs]  # Calculate MAF as minor/major allele count
        maf_average = np.average(mafs)  # Compute the average MAF
        depth_test_maf = maf_average <= min_maf_average  # Check if average MAF meets the threshold
        return depth_test_maf
    else:
        return True  # Skip sites where the major allele frequency is zero

# Function that applies all filtering tests (depth and MAF) to a site.
# If the site passes any test, it will be discarded (return True).
def do_all_tests(depths, allele_counts, max_depths):
    test_min_DP = test_min_depth(depths,)  # Test for minimum depth
    test_max_DP = test_max_depth(depths, max_depths)  # Test for maximum depth
    test_min_mafreq = test_min_maf(allele_counts, max_depths, min_maf_average)  # Test for MAF
    all_tests = True in [test_min_DP, test_max_DP, test_min_mafreq]  # If any test fails, return True
    return all_tests

# Main function that reads sync file, applies filters, and writes the filtered results to output file.
def main(sync_file, sync_outfile, max_depths, min_maf_average):
    with gzip.open(sync_file, 'rt') as data, gzip.open(sync_outfile, 'wb') as file_to_write:
        for sync_line in data:
            allele_counts, depths = get_AC_depths(sync_line, max_depths)  # Extract allele counts and depths

            # If site passes all tests, write it to output file
            if do_all_tests(depths, allele_counts, max_depths) == False:
                file_to_write.write(sync_line.encode())  # Write the line to output file
            else:
                pass  # If site doesn't pass, skip it

# Command-line interface (CLI) setup to run the script with arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-path", type=pathlib.Path)  # Path argument
    parser.add_argument('-sync_file', type=str)  # Sync file input argument
    parser.add_argument("-min_AF", type=float, nargs='?', const=0.05)  # Minimum allele frequency argument (default 0.05)
    parser.add_argument("-min_coverage", type=int, nargs='?', const=10)  # Minimum coverage argument (default 10)
    parser.add_argument("-depth_file", type=str, nargs='?')  # File containing max depths (optional)
    args = parser.parse_args()  # Parse the arguments

    sync_file = args.sync_file

    # Set the minimum depth if provided, else use default
    min_depth = args.min_coverage if args.min_coverage is not None else 10

    # Read max depths from a file if provided, else generate default depths
    if args.depth_file is not None:
        max_depths = pd.read_csv('depths', header=None).iloc[0].tolist()  # Read depths from file
        n_pops = len(max_depths)  # Determine the number of populations
    else:
        n_pops, max_depth = determine_num_pops(sync_file)  # Use default method to determine max depths

    min_maf_average = args.min_AF  # Set minimum MAF average from argument

    # Construct output file name based on input sync file and min MAF value
    sync_outfile = args.sync_file[:-8]+"_DP_MAF" + str(min_maf_average)[1:]+".sync.gz"
    main(sync_file, sync_outfile, max_depths, min_maf_average)  # Run the main processing function

"""
Run this using the following format:

python3 flt_sync_ALT_v5.py \
-sync_file sard_NOZ_VCF.sync.gz \
-min_AF 0.05 \
-min_coverage 10 \
-depth_file depths &

"""
