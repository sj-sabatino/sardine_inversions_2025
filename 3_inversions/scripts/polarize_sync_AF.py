import numpy as np
from random import sample
import gzip
import sys
import numpy as np
from scipy import stats

"""
This script processes sync files containing allele frequency data for populations.
It polarizes allele frequencies based on a specified outgroup and outputs the frequencies
for downstream analysis.

Major Steps:
1. Identify the major allele in the outgroup(s).
2. Polarize allele frequencies for all populations to the major allele of the outgroup.
3. Save the computed allele frequencies to a file.

It outputs a file sardina_polarized_AF.fq that can be used to calculate the AF for
SNPs that are diagnostic for alternative inversion haplotypes. 

"""

# Function to extract chromosome, position, and allele counts from a line in the sync file.
def get_allele_counts(line):
    """
    Parses a sync file line to extract chromosome, position, and allele counts.
    Args:
    - line (str): A line from the sync file.
    Returns:
    - tuple: Chromosome (str), position (str), and allele counts (list of lists).
    """
    line_split = line.strip().split("\t")
    chrm = line_split[0]  # Chromosome
    pos = line_split[1]   # Position
    just_pops = line_split[3:]  # Population data (allele counts)
    allele_counts = [(list(map(int, just_pops[i].split(':')[:4]))) for i in range(0, len(just_pops))]
    return chrm, pos, allele_counts

# Function to determine the major allele across outgroups.
def get_maj_allele_outgroups(allele_counts, outgroup_indicies):
    """
    Finds the major allele in the specified outgroup populations.
    Args:
    - allele_counts (list of lists): Allele counts for all populations.
    - outgroup_indicies (list of int): Indices of the populations to use as outgroups.
    Returns:
    - tuple: Index of the major allele (int) and its corresponding base (str).
    """
    outgroups = [allele_counts[x] for x in outgroup_indicies]  # Select only outgroup populations
    outgroup_maj_indicies = []
    for outgroup in outgroups:
        outgroup_maj_indicies.append(outgroup.index(max(outgroup)))  # Find index of the max allele count
    calc_mode = stats.mode(outgroup_maj_indicies)  # Find the mode of the major alleles
    outgroup_mode = calc_mode.mode[0]
    base_dict = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}  # Mapping allele index to base
    outgroup_base = base_dict[outgroup_mode]
    return outgroup_mode, outgroup_base

# Function to calculate allele frequencies polarized to the outgroup's major allele.
def get_allele_frequencies(line, outgroup_mode):
    """
    Computes allele frequencies relative to the outgroup's major allele.
    Args:
    - line (str): A line from the sync file.
    - outgroup_mode (int): Index of the major allele for polarization.
    Returns:
    - list: Polarized allele frequencies for each population.
    """
    line_split = line.strip().split("\t")
    just_pops = line_split[3:]  # Extract population data
    allele_counts = [(list(map(int, just_pops[i].split(':')[:4]))) for i in range(0, len(just_pops))]
    allele_frequencies = []
    for pop in allele_counts:
        maj_allele_count = int(pop[outgroup_mode])  # Count of the major allele
        del pop[outgroup_mode]  # Remove the major allele from consideration
        other_allele_count = sorted(pop)[-1]  # Largest count among remaining alleles
        freq = maj_allele_count / (maj_allele_count + other_allele_count)  # Compute frequency
        allele_frequencies.append(freq)
    return allele_frequencies

# Main function to process the sync file and save results.
def main(sync_file, sync_outfile):
    """
    Processes a sync file to compute allele frequencies polarized to an outgroup.
    Args:
    - sync_file (str): Path to the input sync file.
    - sync_outfile (str): Path to the output file.
    """
    with open(sync_file, 'r') as data, open(sync_outfile, 'w') as file_to_write:
        for sync_line in data:
            line = sync_line.strip()
            chrom, pos, allele_counts = get_allele_counts(line)
            outgroup_mode, outgroup_base = get_maj_allele_outgroups(allele_counts, outgroup_indicies)
            allele_frequencies = get_allele_frequencies(line, outgroup_mode)
            allele_frequencies = ["{:.5f}".format(x) for x in allele_frequencies]  # Format frequencies
            out_line = "\t".join([chrom, str(pos), outgroup_base] + [str(x) for x in allele_frequencies]) + "\n"
            file_to_write.write(out_line)

# Script entry point
if __name__ == "__main__":
    pop_index = 16 #column position of SMO_J(position 16) in the sync file; AEG_A (positon 1) was also used for chromosome 18 
    outgroup_indicies = [pop_index]
    sync_file = 'sardina_filtered.sync'  # Input sync file
    sync_outfile = 'sardina_polarized_AF.fq'  # Output file
    main(sync_file, sync_outfile)

"""
Example output format:
Chromosome    Position    Major_Allele    Pop1_AF    Pop2_AF    ...
NC_084994.1   9680        G               0.69388    0.26829    ...
"""