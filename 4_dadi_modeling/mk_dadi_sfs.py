#!/usr/bin/env python3
"""
Population genetics script for processing SFS files and SNP data.
Supports a variable number of population files.

Requires:
- A .snps file with format: "scaffold1 1229 A T"
- Multiple SFS files with matching SNPs and no missing data
"""

import sys
import os
from typing import List, Dict


def read_snp_file(path: str, filename: str) -> Dict:
    """Read SNP file and return dictionary of SNP data."""
    snp_dict = {}
    try:
        with open(os.path.join(path, filename), 'r') as infile:
            for line in infile:
                line_split = line.strip().split('\t')
                if len(line_split) < 4:  # Check for proper format
                    continue
                key = f"{line_split[0]}_{line_split[1]}"
                snp_dict[key] = [line_split[2], line_split[3]]
        return snp_dict
    except FileNotFoundError:
        print(f"Error: SNP file {os.path.join(path, filename)} not found.")
        sys.exit(1)


def process_population(pop_data: List[str], maj_allele: str):
    """Process a single population's data."""
    # Check if we have enough columns
    if len(pop_data) < 5:
        return None, None, 0

    try:
        if pop_data[2] == maj_allele:
            maj_count = str(int(float(pop_data[3])))
            min_count = str(int(float(pop_data[4])))
        else:
            maj_count = str(int(float(pop_data[4])))
            min_count = str(int(float(pop_data[3])))
        sum_alleles = int(float(pop_data[3])) + int(float(pop_data[4]))
        return maj_count, min_count, sum_alleles
    except (ValueError, IndexError):
        return None, None, 0


def main():
    # Configuration - can be modified to include any number of populations
    sfs_files = [
        #"Pool1_paired.sfs",  # pop1_AEG_A
        #"Pool2_paired.sfs",  # pop2_AEG_J
        "Pool7_paired.sfs",  # pop3_NPO_A
        #"Pool11_paired.sfs", # pop4_CAN_A
        "Pool16_paired.sfs", # pop5_ION_A
        #"Pool26_paired.sfs"  # pop6_SMO_J
    ]

    # Optional: Define population names
    pop_names = [f"pop{i+1}" for i in range(len(sfs_files))]

    min_count_per_pop = 30
    outfile_name = 'npo_ion_variable.dadi'
    path_to_snpfile = '/home/stephen.sabatino/sardine/popool/vcf/ncbi/snpeff/'
    snp_file = "sardine.snps"

    # Get SNPs
    snp_dict = read_snp_file(path_to_snpfile, snp_file)
    print(f"Loaded {len(snp_dict)} SNPs from file.")

    # Create dynamic header based on number of populations
    pop_header = "\t".join(pop_names)
    header = f"Ingroup\tOutgroup\tAllele1\t{pop_header}\tAllele2\t{pop_header}\tGene\tPosition\n"

    # Initialize counters for reporting
    processed_lines = 0
    filtered_lines = 0
    written_lines = 0

    # Process files
    try:
        # Open all input files
        file_handles = []
        for filename in sfs_files:
            try:
                file_handles.append(open(filename, 'r'))
            except FileNotFoundError:
                print(f"Error: Could not open {filename}")
                # Close already opened files
                for handle in file_handles:
                    handle.close()
                sys.exit(1)

        # Open output file
        with open(outfile_name, 'w') as outfile:
            outfile.write(header)

            # Process files line by line using zip
            for lines in zip(*file_handles):
                processed_lines += 1

                # Parse population data for all populations
                all_pop_data = [line.strip().split('\t') for line in lines]

                # Skip if any line has asterisks in the data
                if any('*' in pop[3] if len(pop) > 3 else False for pop in all_pop_data):
                    filtered_lines += 1
                    continue

                # Skip if any population has insufficient columns for index 7
                if any(len(pop) <= 7 for pop in all_pop_data):
                    filtered_lines += 1
                    continue

                # Skip if any population has empty values at index 7
                if any(pop[7] == '' if len(pop) > 7 else True for pop in all_pop_data):
                    filtered_lines += 1
                    continue

                # Get SNP key from first population and check if it exists
                if len(all_pop_data[0]) < 2:
                    filtered_lines += 1
                    continue

                snp_key = f"{all_pop_data[0][0]}_{all_pop_data[0][1]}"
                if snp_key not in snp_dict:
                    filtered_lines += 1
                    continue

                maj_allele = snp_dict[snp_key][0]
                min_allele = snp_dict[snp_key][1]

                # Process each population
                all_maj_counts = []
                all_min_counts = []
                all_sums = []

                for pop_data in all_pop_data:
                    maj_count, min_count, sum_alleles = process_population(pop_data, maj_allele)
                    if maj_count is None:  # Skip if processing failed
                        break

                    all_maj_counts.append(maj_count)
                    all_min_counts.append(min_count)
                    all_sums.append(sum_alleles)

                # If any population processing failed, skip this SNP
                if len(all_maj_counts) != len(all_pop_data):
                    filtered_lines += 1
                    continue

                # Check if all populations meet minimum count requirement
                if all(sum_alleles > min_count_per_pop for sum_alleles in all_sums):
                    # Join all major counts and minor counts
                    maj_counts_str = "\t".join(all_maj_counts)
                    min_counts_str = "\t".join(all_min_counts)

                    # Write output line
                    line = f"-{maj_allele}-\t---\t{maj_allele}\t{maj_counts_str}\t{min_allele}\t{min_counts_str}\t{snp_key}\t{all_pop_data[0][1]}\n"
                    outfile.write(line)
                    written_lines += 1
                else:
                    filtered_lines += 1

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    finally:
        # Ensure all files are closed
        for handle in file_handles:
            if not handle.closed:
                handle.close()

    # Report statistics
    print(f"Processing complete!")
    print(f"Total lines processed: {processed_lines}")
    print(f"Lines filtered out: {filtered_lines}")
    print(f"Lines written to output: {written_lines}")


if __name__ == "__main__":
    main()
