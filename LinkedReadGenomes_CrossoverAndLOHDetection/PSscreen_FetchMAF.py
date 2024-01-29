#! /usr/bin/env python

### 20211208
### Python script to add minor allelic frequencies for each samples to files of screened phase switched (i.e., putative crossover recombination events)

### Example:
### python ./src/PSscreen_FetchMAF.py ./results/Comparisons/NOV_T504_v_NOV_T501/chr11_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count.bed T504 ./data/Samples/NOV_T504/hapblock_chr11.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf.table_AD_MAF T501 ./data/Samples/NOV_T501/hapblock_chr11.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf.table_AD_MAF

### Import the sys module to access command line arguments
import sys

### Function to read a CSV file and return a list of lists
### Each inner list represents a row, split by the specified split character (default is tab)
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Read the phase switch file specified as the first command line argument
PS=read_csv(sys.argv[1])

### Read the sample names and corresponding Minor Allele Frequency (MAF) files
### These are passed as command line arguments
Sample1=sys.argv[2]
MAF1=read_csv(sys.argv[3])[1:]
Sample2=sys.argv[4]
MAF2=read_csv(sys.argv[5])[1:]

### Determine the output file name based on the input file extension
if '.bed' in sys.argv[1]:
	OutFile=open(sys.argv[1].split('.bed')[0]+'_MAF.bed','w')
else:
	OutFile=open(sys.argv[1]+'_MAF.bed','w')

### Write header line to the output file
OutFile.write(PS[0][0])
for i in range(1,len(PS[0])):
	if PS[0][i]!='INDV':
		OutFile.write('\t'+PS[0][i])
### Include in this header line four new columns, filled with the minor allele frequencies of each sample at the two SNPs between which the putative crossover occurred.
OutFile.write('\t'+Sample1+'_SNP1\t'+Sample1+'_SNP2\t'+Sample2+'_SNP1\t'+Sample2+'_SNP2')

### Iterate through each "phase switch" (or putative crossover)
for p in range(1,len(PS)):
	### If the line has the same length as the first line (this filters out empty rows that occur at the end of some files)
	if len(PS[p])==len(PS[0]):
		### Write the line to the output file
		OutFile.write('\n'+PS[p][0])
		for i in range(1, len(PS[p])):
			OutFile.write('\t'+str(PS[p][i]))
		### Then instantiate a list, creatively called "List"
		### This will hold the positions of the two SNPs between which the phase switch (or putative crossover recombination event) occurred
		List=[]
		### Add the two SNP positions to the list
		List.append(PS[p][1])
		List.append(PS[p][2])
		### Loop through the table of minor allele frequencies for each SNP
		for m in range(0,len(MAF1)):
			### for each SNP position in the list
			for Site in List:
				### again make sure that the line in MAF file isn't an empty line
				if len(MAF1[m])==len(MAF1[0]):
					### Check if the loci in the two files are on the same chromosomal scaffolds AND if their position number is the same
					if MAF1[m][0]==PS[p][0] and MAF1[m][1]==Site:
						### if so, write out the minor allele frequency
						OutFile.write('\t'+str(MAF1[m][3]))
		### Repeat the process for MAF2
		for m in range(0,len(MAF2)):
			for Site in List:
				if len(MAF2[m])==len(MAF2[0]):
					if MAF2[m][0]==PS[p][0] and MAF2[m][1]==Site:
						OutFile.write('\t'+str(MAF2[m][3]))
### Close the output file
OutFile.close()
print('Done :~]')