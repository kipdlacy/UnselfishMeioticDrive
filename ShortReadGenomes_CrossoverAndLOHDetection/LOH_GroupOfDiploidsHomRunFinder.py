#! /usr/bin/env python

# 20210228
# This script is designed to read in a VCF table and for each sample identify all runs of homozygous markers,
# then find the set of unique runs out of those identified in all samples, and then write to an output file
# whether each unique run is present in each sample.

# Usage
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py VCFtableName OutDir

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
# print('SNPhead')
# print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
# print('SNPs[0]')
# print(SNPs[0])
OutDir=sys.argv[2]

# OutFile
# OutFile=open(OutDir+SNPfile.split('/data/SNPtable/')[-1]+"_UniqLOHPresAbs.txt",'w')
OutFile=open(SNPfile+"_UniqLOHPresAbs.txt",'w')
OutFile.write('Chrom\tContig\tLOHstart\tLOHstop\tNumSNPsInLOH\tInfLOHLength\tFirstOnTig\tLastOnTig')
for i in range(6,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# Loop to Identify HomRuns (i.e. LOH) in each sample
print("Identifying HomRuns i.e. LOH in each sample")
LOHlist=[]
HomRun='No'
# For Sample in Table
for i in range(6,len(SNPs[0])):
	print(SNPhead[i])
	LOHlist.append([])
	# For Marker in Table
	for k in range(0,len(SNPs)):
		# If the sample is homozygous at this marker
		if SNPs[k][i][0] == SNPs[k][i][2]:
			# If there is NOT an ongoing run of homozygosity
			if HomRun=='No':
				# Instantiate a new HomRun
				HomList=[SNPs[k][1]]
				HomChrom=SNPs[k][0]
				HomTig=SNPs[k][2]
				# if this marker is on a different contig than the previous marker
				if k!=0 and SNPs[k-1][2] != HomTig:
					FirstOnTig='Yes'
				else:
					FirstOnTig='No'
				HomRun='Yes'
			# Else if there IS an ongoing run of homozygosity
			elif HomRun=='Yes':
				# if current marker is on the same chromosome and conitg as the ongoing LOH run
				if SNPs[k][0]==HomChrom and SNPs[k][2]==HomTig:
					HomList.append(SNPs[k][1])
				# if current marker is on a different chromosome and contig as the ongoing LOH run
				else:
					# If previous marker was on a different contig
					if k!=0 and SNPs[k-1][2]!=SNPs[k][2]:
						LastOnTig='Yes'
					else:
						LastOnTig='No'
					# Append a string describing the finished LOH run to the LOHlist
					LOHlist[-1].append(HomChrom+'-'+HomTig+'-'+HomList[0]+'-'+HomList[-1]+'-'+str(len(HomList))+'-'+str(int(HomList[-1])-int(HomList[0]))+'-'+str(FirstOnTig)+'-'+str(LastOnTig))
					#Instantiate a new HomRun
					HomList=[SNPs[k][1]]
					HomChrom=SNPs[k][0]
					HomTig=SNPs[k][2]
					# if this marker is on a different contig than the previous marker
					if k!=0 and SNPs[k-1][2] != HomTig:
						FirstOnTig='Yes'
					else:
						FirstOnTig='No'
					HomRun='Yes'
		# Else if the sample is heterozygous at this marker
		elif SNPs[k][i][0] != SNPs[k][i][2]:
			# If there IS an ongoing run of homozygosity
			if HomRun=='Yes':
				# If previous marker was on a different contig
				if k!=0 and SNPs[k-1][2]!=SNPs[k][2]:
					LastOnTig='Yes'
				else:
					LastOnTig='No'
				# Append a string describing the finished LOH run to the LOHlist
				LOHlist[-1].append(HomChrom+'-'+HomTig+'-'+HomList[0]+'-'+HomList[-1]+'-'+str(len(HomList))+'-'+str(int(HomList[-1])-int(HomList[0]))+'-'+str(FirstOnTig)+'-'+str(LastOnTig))
				#Deactivate HomRun
				HomRun='No'

# Identify the set of unique HomRuns
Uniques=[]
for Sample in LOHlist:
	for Run in Sample:
		Uniques.append(Run)
Uniques=sorted(set(Uniques))
print('Number of Unique LOH runs: '+str(len(Uniques)))

# Write to the output file
for Run in Uniques:
	OutFile.write('\n'+Run.split('-')[0])
	for i in range(1,len(Run.split('-'))):
		OutFile.write('\t'+Run.split('-')[i])
	for Sample in LOHlist:
		if Run in Sample:
			OutFile.write('\t1')
		else:
			OutFile.write('\t0')

OutFile.close()
print('Done')

