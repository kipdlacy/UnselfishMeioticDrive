#! /usr/bin/env python

# 20210630
# Script to count the number and proportion of sliding windows
# for which one or more genomic features (in this case SNPs)
# are present

### Import Modules
import sys
import numpy

### Defining Functions

# Function to read data from a CSV file
# filename: Path to the CSV file
# l_char: Character to split lines (default: newline character)
# spl_char: Character to split fields within a line (default: tab character)
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments

# Take the first argument from the command line as the input file name
InFile=sys.argv[1]

# Read the data from the input file using the read_csv() function defined above
List=read_csv(InFile)

### Saving the read depth at all sites to a list 
### (note that this is memory intensive: this could be improved in the future)

# Create empty list to store depths
depths=[]
# Iterate through all items in the list (excluding the last two for technical reasons)
for i in List[:-2]:
	# Checking if the length of the list entry is the same as the first
	# This check avoids trying to grab values from entries that are actually from empty lines
	if len(i)==len(List[0]):
		# Append the depth value (converted to an integer) to the list of depths
		depths.append(int(i[2]))

### Calculating Summary Stats

# Converting the depths list to a numpy array for easy statistical calculations
depths=numpy.array(depths)
# calculate various percentiles
quartiles = numpy.percentile(depths, [25, 50, 75])
ninety = numpy.percentile(depths, [5, 95])
ninetyfive = numpy.percentile(depths, [2.5, 97.5])
ninetynine = numpy.percentile(depths, [0.5, 99.5])
# calculate min/max
depths_min, depths_max = depths.min(), depths.max()
# Calculate mean and StdDev
depths_mean=numpy.mean(depths)
depths_std=numpy.std(depths)

### Writing to Output File

# Write the output file name using the input file name
OutFile=open(InFile.split('.depth')[0]+'summary.txt','w')

# Write the header and the computed statistics to the output file
OutFile.write('SampleName\tMin\tQ1\tMedian\tQ3\tMax\t90lo\t90hi\t95lo\t95hi\t99lo\t99hi\tMean\tStdDev\n')
OutFile.write(InFile.split('.depth')[0]+'\t'+str(depths_min)+'\t'+str(quartiles[0])+'\t'+str(quartiles[1])+'\t'+str(quartiles[2])+'\t'+str(depths_max)+'\t'+str(ninety[0])+'\t'+str(ninety[1])+'\t'+str(ninetyfive[0])+'\t'+str(ninetyfive[1])+'\t'+str(ninetynine[0])+'\t'+str(ninetynine[1])+'\t'+str(depths_mean)+'\t'+str(depths_std)+'\n')

# Close the output file
OutFile.close()

