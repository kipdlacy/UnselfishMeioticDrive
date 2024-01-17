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
Count=0
PrevPhase=[0]
for i in range(7,len(SNPs[0])):
	if SNPs[0][i]==SNPs[0][6]:
		PrevPhase.append(0)
	else:
		PrevPhase.append(1)
BlockChrom=SNPs[0][0]
BlockTig=SNPs[0][2]
BlockList=[SNPs[0][1]]
FirstOnTig='Yes'
LastOnTig='No'
for k in range(1,len(SNPs)):
	CurrPhase=[0]
	for i in range(7,len(SNPs[k])):
		if SNPs[k][i]==SNPs[k][6]:
			CurrPhase.append(0)
		else:
			CurrPhase.append(1)
	if SNPs[k][2]!=SNPs[k-1][2]:
		LastOnTig='Yes'
		Count+=1
		OutFile.write('\n'+BlockChrom+'\t'+BlockTig+'\t'+BlockList[0]+'\t'+BlockList[-1]+'\t'+str(len(BlockList))+'\t'+str(int(BlockList[-1])-int(BlockList[0]))+'\t'+FirstOnTig+'\t'+LastOnTig)
		for Samp in PrevPhase:
			OutFile.write('\t'+str(Samp))
		BlockChrom=SNPs[k][0]
		BlockTig=SNPs[k][2]
		BlockList=[SNPs[k][1]]
		FirstOnTig='Yes'
		LastOnTig='No'
	elif CurrPhase!=PrevPhase:
		Count+=1
		OutFile.write('\n'+BlockChrom+'\t'+BlockTig+'\t'+BlockList[0]+'\t'+BlockList[-1]+'\t'+str(len(BlockList))+'\t'+str(int(BlockList[-1])-int(BlockList[0]))+'\t'+FirstOnTig+'\t'+LastOnTig)
		for Samp in PrevPhase:
			OutFile.write('\t'+str(Samp))
		BlockChrom=SNPs[k][0]
		BlockTig=SNPs[k][2]
		BlockList=[SNPs[k][1]]
		FirstOnTig='No'
		LastOnTig='No'
	else:
		BlockList.append(SNPs[k][1])
	PrevPhase=CurrPhase

OutFile.close()
print('Number of Phase Blox: '+str(Count))
print('Done')
