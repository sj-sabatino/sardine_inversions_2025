"""
Steps taken to estimate the allele frequencies for chromosome-scale inversions found
in Sardina pilchardus using  poolseq data from Sabatino et al. 2025.

1. Generate differences in allele frequencies (pwc file) for all SNPs in the filterd sync file for allele
   pairwise comparisions for all populations using Popoolations2
2. Using the PWC file from Popollations, sort "differences in allele frequencies (DAF)" it for the pair of populations with
   that have contrasting allele frequencies for each inverion. Get ths SNPs that are in the top 25th.
3. Generate allele frequencies for all SNPs that are polarized to the population with a high frequency of alternative alleles
   for the inversion of interest.
4. Subset the full SNP allele frequency dataset for the SNP with the highe DAF and may be diagnositc for alternative
   inverstion haplotypes.

"""

# 1 #######################################################################################
# The filtered sync file of allele frequencies was fed into snp-frequency-diff.pl from
# Popoolations2 as follows.
nohup perl snp-frequency-diff.pl \
--input sardina_pilchardus_filtered.sync \
--output-prefix sardina_pilchardus_filtered.allfreq &
###########################################################################################


# 2 #######################################################################################
# Run this for each inversion substiuting the appropriate chromosome and start and stop of the
# putative inverison
# 284 used in the pwc_subset_sort corresponds to the column of pwc for SMO_J population (southern Morocco)
# to which the AF were polarized. SMO_J/CAN_J (column 284 in the pwc file) was used for all inversions except INV_18 which was polarized to the AEG_A/AEG_J
# populations (column 8)
"""
Coordinates for each Inversion used in the small python script just below.

chrm	      start	     stop	       Name
NC_084998.1	0	        32000000	   INV_5
NC_085002.1	1000000	  35108784	   INV_9
NC_085004.1	4000000	  32500000	   INV_11
NC_085005.1	3000000	  34000000	   INV_12
NC_085006.1	16000000	33000000	   INV_13
NC_085008.1	3260000	  29760000	   INV_15
NC_085009.1	3260000	  30000000	   INV_16
NC_085010.1	5500000	  30930000	   INV_17
NC_085011.1	2000000	  4100000	     INV_18
"""

import pandas as pd

#change for each run
chrm = "NC_085001.1"
start = 3260000
stop = 29760000
inversion_name = "INV_5"
pwc_column_position = 284 #corresponds to SMO_J/CAN_J; 8 also used for chrm 18 for AEG_A/AEG_J

pwc=pd.read_csv('sardina_pilchardus_filtered.allfreq_pwc', sep='\t', header=None)
pwc_subset=pwc[(pwc[0]==chrm) & pwc[1]>int(start)) & (pwc[1]<int(stop))]
pwc_subset_sort = pwc_subset.sort_values(pwc_column_position,ascending=False).head(10000).to_csv('snps_to_get_'+inversion_name, sep = '\t', index=False)

#The AF for these ranged from 0 to 1.00, so this gets the top 25th percentile for allele frequency differences
# 284 used in the pwc_subset_sort corresponds to the column of pwc for SMO_J population (southern Morocco)
# to which the AF were polarized.
import pandas as pd
infile = 'snps_to_get_INV5' #change for each run
snps_to_get=pd.read_csv(infile, sep='\t')
snps_to_get_subset75=x[x[pwc_column_position]>=0.75]
pd.concat([snps_to_get_subset75['0'],snps_to_get_subset75['1']], axis=1).to_csv('allele_fq_to_get_'+inversion_name, sep='\t', index=False, header=None)
###########################################################################################

# 3 #######################################################################################
# Make allele frequency file with AF that are polarized to the outgroup. In this case either SMO_J (position 284) for
# all inversions except 18 where AEG_A (8) was used.
polarize_sync_AF.py
###########################################################################################

# 4 #######################################################################################
#Calculate the AF for alternative SNPs in the inversion region.
calc_INV_AF.py
###########################################################################################
