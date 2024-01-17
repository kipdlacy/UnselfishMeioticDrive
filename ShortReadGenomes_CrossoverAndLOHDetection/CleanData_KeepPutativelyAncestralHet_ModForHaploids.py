#! /usr/bin/env python

# 20210227
# Writing script to read in VCF-style table and write only sites that were putatively ancestrally heterozygous in a line to an output file.
# Putative ancestral heterozygosity at any single locus is defined as follows:
# 	Heterozygous in at least one diploid
# 	All heterozygous diploid have the same genotype
# 	Any homozygous diploids bear one of the alles present in the heterozygous genotype

# Usage
# python ./src/CleanData_KeepPutativelyAncestralHet.py VCFtableName

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
OutFile=open(SNPfile+"_PutAncHet",'w')
OutFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

NotFile=open(SNPfile+"_NOTPutAncHet",'w')
NotFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	NotFile.write('\t'+SNPhead[i])

DipInd=[]
HapInd=[]

for i in range(5,len(SNPhead)):
	if 'SM' in SNPhead[i] or 'KL' in SNPhead[i]:
		HapInd.append(i)
		print('Hap: '+SNPhead[i])
	else:
		DipInd.append(i)
		print('Dip: '+SNPhead[i])

for Locus in SNPs:
# Hopefully unneeded troubleshooting
# 	if Locus[0]=='Chr5':
# 		if Locus[1]==3968197 or Locus[1]==3982741:
	if len(Locus)==len(SNPs[0]):
		Hets=[]
		Homs=[]
		HetTru=0
		WriteOut='No'
		for k in DipInd:
			if Locus[k][0] != Locus[k][2]:
				Hets.append(Locus[k])
			else:
				Homs.append(Locus[k])
		if len(Hets) != 0:
			if len(Hets) == 1:
				HetTru+=1
			elif len(Hets) > 1:
				if len(set(Hets))==1:
					HetTru+=1
		if HetTru==1:
			HapPassCount=0
			for k in HapInd:
				if Locus[k][0] == Locus[k][2]:
					if Locus[k][0] == Hets[0][0] or Locus[k][0] == Hets[0][2]:
						HapPassCount+=1
			if HapPassCount==len(HapInd):
				if len(Homs) == 0:
					WriteOut='Yes'
				elif len(Homs) >= 1:
					HomPassCount=0
					for Geno in Homs:
						if Geno[0]==Hets[0][0] or Geno[0]==Hets[0][2]:
							HomPassCount+=1
					if len(Homs)==HomPassCount:
						WriteOut='Yes'
		if WriteOut=='Yes':
			OutFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				OutFile.write('\t'+Locus[i])
		if WriteOut!='Yes':
			NotFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				NotFile.write('\t'+Locus[i])

NotFile.close()
OutFile.close()
print('Done')
