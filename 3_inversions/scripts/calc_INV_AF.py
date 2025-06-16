"""
This script takes a file called sardina_polarized_AF.fq and subsets it for
the SNPs that are diagnostic for alternative inversion haplotypes. 
Needs two files from previous steps: 
1) sardina_polarized_AF.fq
2) allele_fq_to_get
"""

import pandas as pd  # Import pandas

# Define a variable to control the maximum number of SNPs (Single Nucleotide Polymorphisms) to process
num_snps = 5000  # Arbitrarily large number
allele_fq_to_get = "allele_fq_to_get_INV5" # This has to be changed for each inversion. You can make a bash loop etc.

# Initialize an empty list to store SNP identifiers
snp_dict = []

# Open and read the input file containing allele frequency data
with open(allele_fq_to_get, 'r') as infile:
    counter = 0  # Initialize a counter to keep track of processed SNPs
    for i in infile:  
        if counter < num_snps:  # Process only up to the specified number of SNPs
            line_split = i.strip().split("\t")  # Split the line 
            key = line_split[0] + "_" + str(line_split[1])  # Create a unique key using the first two columns
            snp_dict.append(key)  
            counter += 1  
        else:  # Stop processing if the counter exceeds the limit
            break

# Convert the list of SNPs to a set for faster membership checking
snp_dict = set(snp_dict)

# Load a dataset containing SNP data
x = pd.read_csv('sardina_polarized_AF.fq', sep='\t', header=None)

# Create a new column 'key' by combining the first two columns of the dataset
x['key'] = x[0] + "_" + (x[1].astype('str'))  

# Print INV AF
means = list(x[x.key.isin(snp_dict)].mean())  # Filter rows with matching keys and calculate column means
[print(x) for x in means]  # Print the calculated means (column-wise)