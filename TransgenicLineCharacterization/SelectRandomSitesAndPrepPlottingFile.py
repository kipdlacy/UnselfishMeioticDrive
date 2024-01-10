#! /usr/bin/env python

# This script is designed for analyzing DNA sequencing data aligned to two different reference genomes.
# It compares the read depths from a whole-genome sequence aligned to the species' reference genome (WGfile) 
# with those aligned to a transgene insert (TGfile). The known heterozygous variants in the WGfile are used 
# as 'trusted diploid sites' for comparison against normalized read depths in the TGfile.

### Usage
### python ./src/ForR_DepthFileWrangler.py InFileTG InFileWG

### Import Modules
import sys
import random

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
TGfile=sys.argv[1] # Read the file with read depths at each base pair in the transgene insert
TGs=read_csv(TGfile)
WGfile=sys.argv[2] # Read the file with read depths at trusted diploid heterozygous sites
WGs=read_csv(WGfile)

# Creating and naming the output file for combined data
OutFile=open(TGfile+"_TGnWGdepthForR.txt",'w')
# Writing header to the output file
OutFile.write('EndoOrInsert\tChromosome\tPosition\tNormalizedDepth')

# Loop through transgene insert data (TG file)
for Locus in TGs:
	# Writing data to the output file with 'TransgeneInsert' label
	OutFile.write('\nTransgeneInsert\t'+Locus[0]+'\t'+Locus[1]+'\t'+str(Locus[2]))

# Randomly selecting sites from WG file for parallel comparison
# The number of sites selected equals the number of sites in the TG file
for i in sorted(random.sample(range(len(WGs)),len(TGs))):
	# Writing randomly selected data from the WG file with 'TrustedDiploidSites' label
	OutFile.write('\nTrustedDiploidSites\t'+WGs[i][0]+'\t'+WGs[i][1]+'\t'+str(WGs[i][2]))

### Finish
OutFile.close()
print('Done')

