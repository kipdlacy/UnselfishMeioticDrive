#! /usr/bin/env python

### 20220508
### Writing this script to read in a bed file, sum the number of bp included, and divide by the length of the genome composed of chromosomes.

### Import Modules
import sys
import numpy

### Defining a function to read a file and split its content into a list of lists
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Getting the input file name from the command line arguments
InFile=sys.argv[1]
List=read_csv(InFile) # Reading the input file

### Initializing variables to hold the total number of base pairs and SNP count
BpNum=0
SNPnum=0
PhaseBlockLengths = []  # List to hold lengths of each phased block
SNPnumbers = []         # List to hold SNP counts for each block
NumbersOnly = []        # List to hold pairs of block length and SNP count
### Looping through each line in the input file
for Line in List:
	if len(Line) == len(List[1]):
		# Calculating the length of the phased block and updating the total bp count
		BpNum += ((int(Line[2]) - int(Line[1])) + 1)
		# Updating the total SNP count
		SNPnum += int(Line[3])
		# Appending block length and SNP count to respective lists
		PhaseBlockLengths.append((int(Line[2])-int(Line[1]))+1)
		SNPnumbers.append(int(Line[3]))
		NumbersOnly.append([(int(Line[2])-int(Line[1]))+1,int(Line[3])])

### Converting the list to a numpy array and sorting it in descending order based on block length
Array=numpy.array(NumbersOnly)
Array = Array[Array[:, 0].argsort()[::-1]]

### Calculating N50 for SNPs and block lengths
### Initialize variables for the calculation
SNPcount=0
N50snpnum='n/a'
N50blocklength='n/a'
BlockNum=0
### Loop through each block in the sorted array (from largest to smallest)
for Block in Array:
	### Add the number of SNPs in the current block to the total count
	SNPcount+=Block[1]
	### Check if the cumulative SNP count is at least half of the total SNP count
	if SNPcount>=(float(SNPnum)/2):
		### If so, this block is at or crosses the N50 threshold
		### Store the number of SNPs and the length of this block as the N50 values
		N50snpnum=Block[1]
		N50blocklength=Block[0]
		### Break out of the loop as we have found our N50 block
		break
	### Increment the block counter
	BlockNum+=1

### Defining the total size of the genome
GenomeSize=221797443

### Writing the results to an output file
OutFile=open(InFile+"_PropGenomePhased.txt",'w')
OutFile.write('Number of Bp Phased: '+str(BpNum))
OutFile.write('\nNumber of Bp in Reference Genome (Chroms Only): '+str(GenomeSize))
OutFile.write('\nProportion of Genome CO Phased: '+str(float(BpNum)/float(GenomeSize)))
OutFile.write('\nMedian Phase Block Size: '+str(numpy.median(PhaseBlockLengths)))
OutFile.write('\nThere are '+str(SNPnum)+' SNPs in this file, and half of them are found in blocks of size '+str(N50blocklength)+' or longer')

### Printing a summary of results to the console
print(InFile.split('/')[4]+'     '+str(float(BpNum)/float(GenomeSize))+'     '+str(N50blocklength))

### Closing the output file
OutFile.close()
