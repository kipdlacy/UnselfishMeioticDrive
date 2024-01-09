#! /usr/bin/bash

# Here we calculate the depth of DNAseq coverage from aligned files for each sample and obtain summary statistics.

for SAMPLE in Sample1 Sample2 Sample3; do
	# Calculate the depth of coverage using samtools for each position in the BAM file
	# -aa option ensures that all positions are reported, including those with zero depth
	samtools depth -aa ${SAMPLE}_Trimmed_aligned_sorted_merged_dedup.bam >> ${SAMPLE}.depth

	# Run a custom python script to obtain and process depth summary statistics
	python ./scripts/DepthSummaryObtainer.py ./${SAMPLE}.depth
	
	# Remove the depth file after processing to save disk space
	rm ./results/${SAMPLE}.depth
done
