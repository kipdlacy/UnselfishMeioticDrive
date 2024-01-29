#! /usr/bin/env python

### 20211208
### Python script to remove putative phase switches if a previously screened SNP (or any other genomic site with documented issues) lies between the focal SNPs

### Example:
### python ./src/PSscreen_BadSNPsIntervenePS.py ./results/Comparisons/NOV_T504_v_NOV_T501/NOV_T504-v-NOV_T501-chr11_SwitchErrors.diff.switch ./data/Samples/NOV_T504/chr11_AllScreenedSitesUnordered.txt ./data/Samples/NOV_T501/chr11_AllScreenedSitesUnordered.txt

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

### Decision parameter: whether to 'Remove' phase switches with dirty sites or 'Count' them
RemoveOrCount=sys.argv[5]

### If 'Remove' is specified, this block executes
if RemoveOrCount=='Remove':
	### Determine the output file name and open it for writing
	if '.bed' in sys.argv[1]:
		OutFile=open(sys.argv[1].split('.bed')[0]+'_No'+sys.argv[4]+'Inter.bed','w')
	else:
		OutFile=open(sys.argv[1]+'_No'+sys.argv[4]+'Inter.bed','w')
	### Write the header to the output file
	OutFile.write(PS[0][0])
	for i in range(1,len(PS[0])):
		if PS[0][i]!='INDV':
			OutFile.write('\t'+PS[0][i])
	### Process each phase switch (putative crossover recombination event)
	for p in range(1,len(PS)):
		### Check if the line is the same length as the first one (this is to avoid processing empty lines which sometimes occur at the end of files)
		if len(PS[p])==len(PS[0]):
			### Instantiate "WriteOut" variable which is used to hold the information of whether a line has passed screening or not
			WriteOut='Yes'
			### Check each "Dirty" site in the file for Sample 1
			for k in Dirty1:
				### Again check to be sure the line isn't empty
				if len(k)==len(Dirty1[0]):
					### If the SNPs are on the same chromosome, 
					### and if the position of the "dirty" site is greater than the position of the 
					### first SNP in the crossover event but less than the position of the second SNP (i.e., is between or "intervenes" these two)
					if k[0]==PS[p][0] and int(k[1]) > int(PS[p][1]) and int(k[1]) < int(PS[p][2]):
						### Don't write this line to the output file
						WriteOut='No'
			### Same FOR each "Dirty" site in the file for Sample 2
			for k in Dirty2:
				if len(k)==len(Dirty2[0]):
					if k[0]==PS[p][0] and int(k[1]) > int(PS[p][1]) and int(k[1]) < int(PS[p][2]):
						WriteOut='No'
			### Write to the output file if it hasn't satisfied these criteria
			if WriteOut=='Yes':
				OutFile.write('\n'+PS[p][0])
				for j in range(1, len(PS[p])):
					if PS[0][j]!='INDV':
						OutFile.write('\t'+PS[p][j])
	OutFile.close()
### If 'Count' is specified, this block executes
elif RemoveOrCount=='Count':
	### Determine the output file name and open it for writing
	if '.bed' in sys.argv[1]:
		OutFile=open(sys.argv[1].split('.bed')[0]+'_'+sys.argv[4]+'InterCount.bed','w')
	else:
		OutFile=open(sys.argv[1]+'_'+sys.argv[4]+'InterCount.bed','w')
	### Write the header to the output file
	OutFile.write(PS[0][0])
	for i in range(1,len(PS[0])):
		if PS[0][i]!='INDV':
			OutFile.write('\t'+PS[0][i])
	### Add a column to hold the counts
	OutFile.write('\t'+sys.argv[4]+'InterveneCounts')
	### Process each phase switch (putative crossover recombination event)
	for p in range(1,len(PS)):
		### Check if the line is the same length as the first one (this is to avoid processing empty lines which sometimes occur at the end of files)
		if len(PS[p])==len(PS[0]):
			### Instantiate counter
			Count=0
			### Check each "Dirty" site in the file for Sample 1
			for k in Dirty1:
				### Again, avoid empty lines
				if len(k)==len(Dirty1[0]):
					### If the SNPs are on the same chromosome, 
					### and if the position of the "dirty" site is greater than the position of the 
					### first SNP in the crossover event but less than the position of the second SNP (i.e., is between or "intervenes" these two)
					if k[0]==PS[p][0] and int(k[1]) > int(PS[p][1]) and int(k[1]) < int(PS[p][2]):
						### Count it
						Count+=1
			### Again for Sample 2
			for k in Dirty2:
				if len(k)==len(Dirty2[0]):
					if k[0]==PS[p][0] and int(k[1]) > int(PS[p][1]) and int(k[1]) < int(PS[p][2]):
						Count+=1
			### Write the line plus the count to the output file
			OutFile.write('\n'+PS[p][0])
			for j in range(1, len(PS[p])):
				if PS[0][j]!='INDV':
					OutFile.write('\t'+PS[p][j])
			OutFile.write('\t'+str(Count))
	OutFile.close()

