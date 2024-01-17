#! /usr/bin/env python

# 20210226
# Copied from ~/scripts/crap/Planning/20210217*/src/ReadDepth_AppendContigNamesToSNPs_20210219.py

# 20210219
# This script is designed to read in a file of snps and a file of contig breaks and ask whether any snps fall within 1000bp of the beginning or end of a contig, writing such snps to an output file.

# USage
# python ./src/ReadDepth_AppendContigNamesToSNPs_20210219.py ./data/SNPs_NoSVs_20210219.txt ../../Papers/WingedMutant/data/SlidingWindows/Input/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
print('SNPhead:')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]:')
print(SNPs[0])
ContigFile=sys.argv[2]
ContigHead=read_csv(ContigFile)[0]
print('ContigHead:')
print(ContigHead)
Contigs=read_csv(ContigFile)[1:]
print('Contigs[0]:')
print(Contigs[0])

# Creating Output file
OutFile=open(SNPfile+"_Contigs",'w')
OutFile.write(str(SNPhead[0])+'\t'+str(SNPhead[1])+'\t'+str(ContigHead[3]))
for i in range(2,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# Double for loop
NumSNPsNearBreaks=0
TotalNumSNPs=0
for SNP in SNPs:
	for Wind in Contigs:
		#integerizing the site
		Site=int(SNP[1])
		#now checking if its the same locus
		if SNP[0]==Wind[0]:
			if Site >= int(Wind[1]) and Site <= int(Wind[2]):
				OutFile.write('\n'+SNP[0]+'\t'+str(Site)+'\t'+str(Wind[3]))
				for i in range(2,len(SNP)):
					OutFile.write('\t'+SNP[i])

OutFile.close()
print('done :)!')
