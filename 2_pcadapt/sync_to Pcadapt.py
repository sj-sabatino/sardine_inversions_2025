#!/usr/bin/env python3

"""
This script takes a sync file and creates the input data necessary for PCAdapt.
It implements a sliding window approach with a given number of SNPs per window and step size
"""

import gzip
from collections import defaultdict
import pickle
import numpy as np
import argparse

# Class to represent a SNP (Single Nucleotide Polymorphism)
class SNP(object):
    def __init__(self, pos, freq, major, minor):
        """
        Initialize a SNP object with position, frequencies, major and minor alleles.
        :param pos: Position of the SNP (formatted as a string)
        :param freq: List of frequencies for the minor allele in different samples
        :param major: The major allele (reference allele)
        :param minor: The minor allele (variant allele)
        """
        self.pos = pos
        self.freq = freq
        self.major = major
        self.minor = minor

    def __repr__(self):
        """
        String representation for the SNP object.
        :return: Formatted string describing the SNP details
        """
        return "Position: {}\nMajor allele: {}\nMinor allele: {}\nMinor Allele frequencies: {}\n".format(self.pos, self.major, self.minor, self.freq)


# Function to parse the BED file, returning the scaffolds and SNP positions to consider
def bed_snps(bed_path):
    scaffolds = set()
    used_snps = set()

    # Read the bed file
    with gzip.open(bed_path, 'rt') as bed:
        count = 0
        for line in bed:
            if count == 0:
                count += 1
                continue  # Skip the header line
            data = line.strip().replace("\t", "_")  # Replace tabs with underscores to create a unique identifier
            scaffolds.add(line.strip().split("\t")[0])  # Add scaffold names
            used_snps.add(data)  # Add SNP positions

    # Sort scaffolds numerically (assuming scaffold names are prefixed with "scaffold" followed by a number)
    scaffolds = sorted(list(scaffolds), key=lambda x: int(x.replace("scaffold", "")))
    return scaffolds, used_snps


# Function to parse a VCF (Variant Call Format) file
def parse_vcf(vcf_path, used_snps):
    snps = defaultdict(SNP)
    # Read the VCF file
    with gzip.open(vcf_path, 'rt') as vcf:
        header = True
        for line in vcf:
            if header:
                header = False
                # Get the number of samples from the VCF header (after the first 8 columns)
                n_samples = len(line.strip().split("\t")[9::])
            else:
                data = line.strip().split("\t")  # Split VCF data by tabs
                pos = "_".join(data[0:2])  # Combine chromosome and position as a unique key
                if pos not in used_snps:
                    continue  # Skip SNPs that are not in the bed file
                where = data[8].split(":")  # Extract allele information from INFO column
                allele1_pos = where.index("RO")
                allele1_base = data[3]  # Reference allele

                allele2_pos = where.index("AO")
                allele2_base = data[4]  # Alternative allele

                fpop_features = data[9].split(":")
                if len(fpop_features) == 1:
                    fpop_features = data[10].split(":")
                    if len(fpop_features) == 1:
                        continue  # Skip SNP if no allele information available

                # Determine the major and minor alleles based on counts
                if int(fpop_features[allele1_pos]) >= int(fpop_features[allele2_pos]):
                    major_pos = allele1_pos
                    major_base = allele1_base
                    minor_pos = allele2_pos
                    minor_base = allele2_base
                else:
                    major_pos = allele2_pos
                    major_base = allele2_base
                    minor_pos = allele1_pos
                    minor_base = allele1_base

                # Calculate the frequency of the minor allele for each sample
                frequencies = []
                for i in range(9, 9+n_samples):
                    features = data[i].split(":")
                    if len(features) == 1:  # Handle missing data
                        frequencies += [0]  # Assign zero for missing data
                        continue
                    major = int(features[major_pos])
                    minor = int(features[minor_pos])
                    if minor == 0:
                        frequencies += [0]  # Assign zero if no minor allele is present
                        continue
                    freq = minor / (major + minor)  # Calculate frequency of minor allele
                    frequencies += [freq]

                # Create SNP object and store in dictionary
                c_SNP = SNP(pos, frequencies, major_base, minor_base)
                snps[pos] = c_SNP
    return snps


# Function to parse a SYNC file
def parse_sync(sync_path, used_snps):
    snps = defaultdict(SNP)
    base_pos = {0: 'A', 1: 'T', 2: 'C', 3: 'G', 4: 'N', 5: 'del'}  # Mapping numeric positions to base pairs
    # Read the SYNC file
    with gzip.open(sync_path, 'rt') as sync:
        for line in sync:
            scaff, bpos, inrefall, *individuals = line.strip().split("\t")
            individuals = [list(map(int, x.split(":")[0:4])) for x in individuals]  # Convert to list of integers
            pos = "{}_{}".format(scaff, bpos)  # Create unique SNP identifier
            if pos not in used_snps:
                continue  # Skip SNPs not in the bed file
            # Calculate the counts for each allele
            counts = [sum(x) for x in zip(*individuals)]
            major_allele = sorted(counts, reverse=True)
            major_pos = counts.index(major_allele[0])
            counts[major_pos] = 0  # Reset major allele position to avoid selecting it again
            minor_pos = counts.index(major_allele[1])
            major_base = base_pos[major_pos]  # Map position to nucleotide base
            minor_base = base_pos[minor_pos]

            # Calculate the frequency of the minor allele for each sample
            frequencies = []
            for i in individuals:
                major = int(i[major_pos])
                minor = int(i[minor_pos])
                if minor == 0:
                    frequencies += [0]  # Assign zero if no minor allele is present
                    continue
                freq = minor / (major + minor)  # Calculate frequency of minor allele
                frequencies += [freq]

            # Create SNP object and store in dictionary
            c_SNP = SNP(pos, frequencies, major_base, minor_base)
            snps[pos] = c_SNP
    return snps


# Function to calculate average frequencies across sliding windows
def average_across_windows(snps, window_size, step_size, n_scaffolds=1000):
    scaffold_snps = defaultdict(list)
    scaffolds = set()
    results = list()
    header = list()
    
    # Group SNPs by scaffolds
    for pos, snp in snps.items():
        scaffold = "_".join(pos.split("_")[0:2])
        scaffolds.add(scaffold)
        scaffold_snps[scaffold] += [snp]

    # Sort scaffolds numerically
    scaffolds = sorted(list(scaffolds), key=lambda x: int(x.split("_")[1]))

    # Process each scaffold and create sliding windows
    for i, scaff in enumerate(scaffolds):
        if i == n_scaffolds:
            break
        snps_to_use = scaffold_snps[scaff]
        i = 0
        while i < len(snps_to_use) - window_size + 1:
            start_pos = snps_to_use[i].pos
            stop_pos = snps_to_use[i + window_size - 1].pos
            snps_freqs = np.array([x.freq for x in snps_to_use[i:i + window_size]])
            means = np.mean(snps_freqs, axis=0, dtype=np.float64)  # Compute mean across the window
            results += [means]
            header += ["{}:{}".format(start_pos, stop_pos)]  # Window start and end positions
            i += step_size  # Move the window by the step size

    return header, results


# Function to process SNPs without using sliding windows
def no_windows(snps, n_scaffolds=1000):
    scaffold_snps = defaultdict(list)
    scaffolds = set()
    results = list()
    header = list()

    # Group SNPs by scaffolds
    for pos, snp in snps.items():
        scaffold = "_".join(pos.split("_")[0:2])
        scaffolds.add(scaffold)
        scaffold_snps[scaffold] += [snp]

    # Sort scaffolds numerically
    scaffolds = sorted(list(scaffolds), key=lambda x: int(x.split("_")[1]))

    # Process each scaffold without windows
    for i, scaff in enumerate(scaffolds):
        if i == n_scaffolds:
            break
        snps_to_use = scaffold_snps[scaff]
        i = 0
        for snp in snps_to_use:
            results += [snp.freq]
            header += ["{}".format(snp.pos)]  # SNP positions

    return header, results


# Function to compare SNP data between two datasets (sync vs vcf)
def pair_snps(snps1, snps2):
    for sync, vcf in zip(snps1.values(), snps2.values()):
        print(sync.pos, vcf.pos)
        print(sync.major, vcf.major)
        print(sync.minor, vcf.minor)
        d = sync.freq
        d = [d[0], d[1], d[3], d[4], d[2]]  # Reorder frequencies
        print(d)
        print(vcf.freq)
        print()
        input()  # Pause to inspect results


# Function to write the sliding window results to an output file
def write_output(windows, header, out_path):
    with open(out_path, 'w') as out:
        n_samples = len(windows[0])
        out.write("Sample\t{}\n".format("\t".join(header)))
        for i in range(n_samples):
            out.write("{}\t{}\n".format(i, "\t".join([str(x[i]) for x in windows])))


# Function to write only the frequencies to an output file (without sample labels)
def write_output_only_freqs(windows, header, out_path):
    with open(out_path, 'w') as out:
        n_samples = len(windows[0])
        out.write("{}\n".format("\t".join(header)))
        for i in range(n_samples):
            out.write("{}\n".format("\t".join([str(x[i]) for x in windows])))


# Main execution starts here
if __name__ == '__main__':
    '''
    usage: makepcaadaptinput.py [-h] -b bed_file -sf sync_file -vf vcf_file -o output -w window_size -ss step_size

    Build PCAdapt files from vcf or sync inputs (a list of snps to consider needs to be provided in a separate bed file)


    Required Arguments:
    -b bed_file         Path to a bed file that contains the list of snps to consider, it expects a header in the form of CHR\tpos
    -sf synf_file       Path to a sync file to build the PCAdapt file from
    -vf vcf_file        Path to a vcf file to build the PCAdapt file from
    -o output           Path to the output folder
    -w window_size      Size (in number of snps) of the window # default = 20
    -ss step_size       Step size (in number of snps) for the sliding window
    -n number_scaffolds Number of scaffolds to use

    '''

    # Argument parser to handle input parameters
    parser = argparse.ArgumentParser(description="Build PCAdapt file from sync or vcf")
    parser.add_argument('-b', metavar='bed_file', type=str, required=True,
                        help='A path to a bed file that contains the list of snps to consider, it expects a header in the form of CHR\tpos')
    parser.add_argument('-sf', metavar='synf_file', type=str, required=False, default=None,
                        help='Path to a sync file to build the PCAdapt file from')
    parser.add_argument('-vf', metavar='vcf_file', type=str, required=False, default=None,
                        help='Path to a vcf file to build the PCAdapt file from')
    parser.add_argument('-o', metavar='output', type=str, required=True,
                        help='Path to the output folder')
    parser.add_argument('-w', metavar="window_size", type=int, default=20,
                        help='Size (in number of snps) of the window # default = 20')
    parser.add_argument('-ss', metavar="step_size", type=int, default=10,
                        help='Step size (in number of snps) for the sliding window # default = 10')
    parser.add_argument('-n', metavar='number_of_scaffolds', type=int, default=1000,
                        help='Number of scaffolds to include, they will be ordered by scaffold name')

    # Parse arguments
    args = parser.parse_args()

    # Ensure at least one input file (sync or vcf) is provided
    assert args.sf is not None or args.vf is not None, "Please provide either a vcf or a sync file"

    # Parse the SNPs from the bed file
    scaffolds, used_snps = bed_snps(args.b)

    # If a sync file is provided, process it
    if args.sf is not None:
        snpssync = parse_sync(args.sf, used_snps)
        sheader, sync_windows = average_across_windows(snpssync, args.w, args.ss, args.n)
        write_output_only_freqs(sync_windows, sheader, "{}sync.pcadapt".format(
            args.o if args.o[-1] == "/" else "{}/".format(args.o)))

    # If a VCF file is provided, process it
    if args.vf is not None:
        snpsvcf = parse_vcf(args.vf, used_snps)
        vheader, vcf_windows = average_across_windows(snpsvcf, args.w, args.ss, args.n)
        write_output_only_freqs(vcf_windows, vheader, "{}vcf.pcadapt".format(
            args.o if args.o[-1] == "/" else "{}/".format(args.o)))

    # If no window is used, process SNPs directly
    if args.sf is not None and args.w == 1:
        print("using no window")
        snpssync = parse_sync(args.sf, used_snps)
        sheader, sync_snps = no_windows(snpssync, args.n)
        write_output_only_freqs(sync_snps, sheader, "{}sync.pcadapt".format(
            args.o if args.o[-1] == "/" else "{}/".format(args.o)))
    
    if args.vf is not None and args.w == 1:
        snpsvcf = parse_vcf(args.vf, used_snps)
        sheader, vcf_snps = no_windows(snpsvcf, args.n)
        write_output_only_freqs(vcf_snps, sheader, "{}sync.pcadapt".format(
            args.o if args.o[-1] == "/" else "{}/".format(args.o)))