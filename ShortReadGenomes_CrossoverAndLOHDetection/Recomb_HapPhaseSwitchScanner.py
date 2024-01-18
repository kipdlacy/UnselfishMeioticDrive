#! /usr/bin/env python

# 20210302
# This script is written to read in a VCF-style table of variants in haploid samples and document evident recombination events among samples. 

# Usage
# python ./src/Recomb_HapPhaseSwitchScanner.py VCFtable

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

# OutFile
OutFile=open(SNPfile+"_PhaseSwitches.txt",'w')
OutFile.write('Chrom\tContig\tBlockStart\tBlockStop\tNumSNPsInBlock\tInfBlockLength\tFirstOnTig\tLastOnTig')
for i in range(6,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# Loop
Count=0 # counting the number of phase blocks
PrevPhase=[0] # instantiating a list to hold the relative phasing of samples
# Index 6 is the first sample in this table, and the rest are in indices 7 through the end
# so this for loop starts with index 7, and for that Sample and each subsequent one
for i in range(7,len(SNPs[0])):
	# compares the genotype to that of the sample in index 6
	if SNPs[0][i]==SNPs[0][6]:
		# if they are the same, the relative phase is given as 0
		PrevPhase.append(0)
	else:
		# if different, relative phase is given as 1
		PrevPhase.append(1)
BlockChrom=SNPs[0][0] # instantiate a variable for the chromosome the current phase block is on
BlockTig=SNPs[0][2] # instantiate a variable for the contig the current phase block is on
BlockList=[SNPs[0][1]] # create a list containing the names of the numerical position of each SNP in the phase block
FirstOnTig='Yes' # the first SNP is the first on the contig
LastOnTig='No' # and is not the last
# now, startingwith the second SNP, and looping across all of them
for k in range(1,len(SNPs)):
	# instantiate a list of the current phase. The sample in the 6th index always has relative phase = 0, because it's genotype is the same as it's own genotype
	CurrPhase=[0]
	# Index 6 is the first sample in this table, and the rest are in indices 7 through the end
	# so this for loop starts with index 7, and for that Sample and each subsequent one
	for i in range(7,len(SNPs[k])):
		# compares the genotype to that of the sample in index 6
		if SNPs[k][i]==SNPs[k][6]:
			# if they are the same, the relative phase is given as 0
			CurrPhase.append(0)
		else:
			# if different relative phase is given as 1
			CurrPhase.append(1)
	# if the contig that the current SNP is on is different than the one the previous SNP was on
	if SNPs[k][2]!=SNPs[k-1][2]:
		# then the last SNP was the last on it's contig
		LastOnTig='Yes'
		# the number of phase blocks is increased
		Count+=1
		# write the information for the block that just was terminated to the output file
		OutFile.write('\n'+BlockChrom+'\t'+BlockTig+'\t'+BlockList[0]+'\t'+BlockList[-1]+'\t'+str(len(BlockList))+'\t'+str(int(BlockList[-1])-int(BlockList[0]))+'\t'+FirstOnTig+'\t'+LastOnTig)
		# then write the relative phasing for that phase block
		for SampRelPhase in PrevPhase:
			OutFile.write('\t'+str(SampRelPhase))
		BlockChrom=SNPs[k][0] # set the new 'BlockChrom' with the info for the current SNP
		BlockTig=SNPs[k][2] # set the new 'BlockTig' with the info for the current SNP
		BlockList=[SNPs[k][1]] # instantiate the new 'BlockList' with the numerical position of the current SNP
		FirstOnTig='Yes' # since the previous tig just ended this is now the first on the contig
		LastOnTig='No' # reinstantiate LastOnTig with a value of 'no'
	# if the relative phasing of the samples has changed
	elif CurrPhase!=PrevPhase:
		# add to the number of phase blocks
		Count+=1
		# write to the outfile
		OutFile.write('\n'+BlockChrom+'\t'+BlockTig+'\t'+BlockList[0]+'\t'+BlockList[-1]+'\t'+str(len(BlockList))+'\t'+str(int(BlockList[-1])-int(BlockList[0]))+'\t'+FirstOnTig+'\t'+LastOnTig)
		for SampRelPhase in PrevPhase:
			OutFile.write('\t'+str(SampRelPhase))
		BlockChrom=SNPs[k][0] # set the new 'BlockChrom' with the info for the current SNP
		BlockTig=SNPs[k][2] # set the new 'BlockTig' with the info for the current SNP
		BlockList=[SNPs[k][1]] # instantiate the new 'BlockList' with the numerical position of the current SNP
		FirstOnTig='No' # since the previous block was not the last on it's contig, this can't be the first on it's contig
		LastOnTig='No' # reinstantiate LastOnTig with a value of 'no'
	# if the relative phasing of the current and previous SNPs are the same, and both are on the same contig
	else:
		# add the numerical position of the current SNP to the list
		BlockList.append(SNPs[k][1])
	# change the current phase to the previous phase before starting the loop to the next SNP
	PrevPhase=CurrPhase

OutFile.close()
print('Number of Phase Blox: '+str(Count))
print('Done')
