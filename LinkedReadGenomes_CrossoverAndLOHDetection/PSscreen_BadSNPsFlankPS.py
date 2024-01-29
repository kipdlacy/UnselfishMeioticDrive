#! /usr/bin/env python

### 20211208
### Python script to remove putative phase switches if a previously screened SNP lies between either focal SNP and the
### SNP two upstream or two downstream.
### This script takes as input 
### 1) a list of phase switches (putative crossover recombination events)
### 2 & 3) Lists of "Dirty" sites for each sample which are associated with genome assembly errors and indicate that phasing cannot be trusted for that region
### 4 & 5) The VCF tables for each sample. These are used to identify the locations of the upstream and downstream SNPs
### 6) The user-input name of the type of site that will be screened against. This is just used for naming the output file.
### 7) The user-input number of SNPs within which the "dirty sites" are not allowed. This is used in the screening process and in naming the file.
### 8) An option of whether to Remove PS from the file, or to count that number of PS SNPs that are screened

### Example:
### python ./src/PSscreen_BadSNPsFlankPS.py ./results/Comparisons/NOV_T504_v_NOV_T501/NOV_T504-v-NOV_T501-chr11_SwitchErrors.diff.switch_NoScreenedSNPsIntervene.bed ./data/Samples/NOV_T504/chr11_AllScreenedSitesUnordered.txt ./data/Samples/NOV_T501/chr11_AllScreenedSitesUnordered.txt ./data/Samples/NOV_T504/hapblock_chr11.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf.table ./data/Samples/NOV_T501/hapblock_chr11.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf.table

### Importing the sys module to use command line arguments
import sys

### Function to read a CSV file. Returns a list of lists, where each sublist is a row in the file.
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Reading the phase switch (PS) data from the first command line argument
PS=read_csv(sys.argv[1])

### Reading the 'dirty' sites data from the next two command line arguments. 
### These are problematic genomic sites that need to be screened out.
### [1:] means that I am skipping the header
Dirty1=read_csv(sys.argv[2])[1:]
Dirty2=read_csv(sys.argv[3])[1:]

### Reading the VCF-style tables of variants for each sample from the next two command line arguments. 
### These are used to identify the genomic positions of each upstream and downstream SNP
### [1:] means that I am skipping the header
VCF1=read_csv(sys.argv[4])[1:]
VCF2=read_csv(sys.argv[5])[1:]

### The number of SNPs upstream or downstream that constitute what "flanking" is
Num=int(sys.argv[7])

### Option of whether to Remove SNPs from the output file, or count
RemoveOrCount=sys.argv[8]

### If 'Remove' is specified, this block executes
if RemoveOrCount=='Remove':
	### Determine the output file name and open it for writing
	if '.bed' in sys.argv[1]:
		OutFile=open(sys.argv[1].split('.bed')[0]+'_No'+sys.argv[6]+'WN_'+sys.argv[7]+'.bed','w')
	else:
		OutFile=open(sys.argv[1]+'_No'+sys.argv[6]+'WN_'+sys.argv[7]+'.bed','w')
	### Write the header line to the output file.
	OutFile.write(PS[0][0])
	for i in range(1,len(PS[0])):
		if PS[0][i]!='INDV':
			OutFile.write('\t'+PS[0][i])
	### Process each phase switch.
	for p in range(1,len(PS)):
		### Check if the line is the same length as the first one (this is to avoid processing empty lines which sometimes occur at the end of files)
		if len(PS[p])==len(PS[0]):
			### Instantiate "WriteOut" variable which is used to hold the information of whether a line has passed screening or not
			WriteOut='Yes'
			### Check for dirty sites around the phase switch in the first sample.
			### First identify the upstream and downstream SNPs in the first sample by looping through the VCF-style genotype table
			for m in range(0,len(VCF1)):
				### again check to make sure you're not checking empty lines
				if len(VCF1[m])==len(VCF1[0]):
					### if the chromosome name is the same and the SNP position is the same as that of the first SNP in the putative crossover event
					if VCF1[m][0]==PS[p][0] and VCF1[m][1]==PS[p][1]:
						### Then loop through the "dirty" file for that sample
						for k in Dirty1:
							### again check to make sure you're not checking empty lines
							if len(k)==len(Dirty1[0]):
								### if the chromosome's the same, and the bad SNP lies between the upstream SNP and the sample
								if k[0]==PS[p][0] and int(k[1]) > int(VCF1[m-Num][1]) and int(k[1]) < int(VCF1[m][1]):
									WriteOut='No'
					### Same but for the second SNP and the downstream SNP
					elif VCF1[m][0]==PS[p][0] and VCF1[m][1]==PS[p][2]:
						for k in Dirty1:
							if len(k)==len(Dirty1[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF1[m][1]) and int(k[1]) < int(VCF1[m+Num][1]):
									WriteOut='No'
			### Repeat for the second sample
			for m in range(0,len(VCF2)):
				if len(VCF2[m])==len(VCF2[0]):
					if VCF2[m][0]==PS[p][0] and VCF2[m][1]==PS[p][1]:
						for k in Dirty2:
							if len(k)==len(Dirty2[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF2[m-Num][1]) and int(k[1]) < int(VCF2[m][1]):
									WriteOut='No'
					elif VCF2[m][0]==PS[p][0] and VCF2[m][1]==PS[p][2]:
						for k in Dirty2:
							if len(k)==len(Dirty2[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF2[m][1]) and int(k[1]) < int(VCF2[m+Num][1]):
									WriteOut='No'
			### Write to the output file
			if WriteOut=='Yes':
				OutFile.write('\n'+PS[p][0])
				for j in range(1, len(PS[p])):
					if PS[0][j]!='INDV':
						OutFile.write('\t'+PS[p][j])
### If 'Remove' is specified, this block executes
### THe logic is the same as above, but we write out all putative crossovers, just with counting the number of matches
elif RemoveOrCount=='Count':
	if '.bed' in sys.argv[1]:
		OutFile=open(sys.argv[1].split('.bed')[0]+'_'+sys.argv[6]+'WN_'+sys.argv[7]+'_Count.bed','w')
	else:
		OutFile=open(sys.argv[1]+'_'+sys.argv[6]+'WN_'+sys.argv[7]+'_Count.bed','w')
	OutFile.write(PS[0][0])
	for i in range(1,len(PS[0])):
		if PS[0][i]!='INDV':
			OutFile.write('\t'+PS[0][i])
	OutFile.write('\tWithin'+sys.argv[7]+'SNPsOfPSCount')
	for p in range(1,len(PS)):
		if len(PS[p])==len(PS[0]):
			Count=0
			for m in range(0,len(VCF1)):
				if len(VCF1[m])==len(VCF1[0]):
					if VCF1[m][0]==PS[p][0] and VCF1[m][1]==PS[p][1]:
						for k in Dirty1:
							if len(k)==len(Dirty1[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF1[m-Num][1]) and int(k[1]) < int(VCF1[m][1]):
									Count+=1
					elif VCF1[m][0]==PS[p][0] and VCF1[m][1]==PS[p][2]:
						for k in Dirty1:
							if len(k)==len(Dirty1[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF1[m][1]) and int(k[1]) < int(VCF1[m+Num][1]):
									Count+=1
			for m in range(0,len(VCF2)):
				if len(VCF2[m])==len(VCF2[0]):
					if VCF2[m][0]==PS[p][0] and VCF2[m][1]==PS[p][1]:
						for k in Dirty2:
							if len(k)==len(Dirty2[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF2[m-Num][1]) and int(k[1]) < int(VCF2[m][1]):
									Count+=1
					elif VCF2[m][0]==PS[p][0] and VCF2[m][1]==PS[p][2]:
						for k in Dirty2:
							if len(k)==len(Dirty2[0]):
								if k[0]==PS[p][0] and int(k[1]) > int(VCF2[m][1]) and int(k[1]) < int(VCF2[m+Num][1]):
									Count+=1
			OutFile.write('\n'+PS[p][0])
			for j in range(1, len(PS[p])):
				if PS[0][j]!='INDV':
					OutFile.write('\t'+PS[p][j])
			OutFile.write('\t'+str(Count))

### Close the output file
OutFile.close()
