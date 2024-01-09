#! /usr/bin/bash

# This file contains the commands I used to process and align standard short read whole genome sequencing data to the O. biroi reference genome
# For simplicity, I don't include the names of all samples used in this stud, and instead loop across hypothetical samples.

# Initialize variable to count samples, used for read group labeling during alignment
SampleCount=1
# Loop through each sample
for SAMPLE in Sample1 Sample2 Sample3; do

	# Trim the raw data using Trimmomatic
	java -jar /Path/To/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
		${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz \
		-baseout ${SAMPLE}_Trimmomatic_Out_ \
		ILLUMINACLIP:/Path/To/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:10:10 \
		HEADCROP:16 SLIDINGWINDOW:4:15 MINLEN:50

	# Align the trimmed reads to the reference genome
	bwa mem -M -R "@RG\tID:${SAMPLE}_group\tSM:${SAMPLE}\tPL:illumina\tLB:lib${SampleCount}\tPU:unit${SampleCount}" \
		/Path/To/Reference/Genome/Obir.assembly.v5.4.fasta \
		${SAMPLE}_Trimmomatic_Out__1P ${SAMPLE}_Trimmomatic_Out__2P \
		> ${SAMPLE}_L1_Trimmed_P_bwa_aligned.sam
	# Increase the SampleCount for the next sample's read group labeling
	SampleCount=$((SampleCount + 1))

	# Sort the aligned reads using
	java -jar /Path/To/picard/picard.jar SortSam \
		INPUT=${SAMPLE}_L1_Trimmed_P_bwa_aligned.sam \
		OUTPUT=${SAMPLE}_Trimmed_P_bwa_aligned_sorted.bam \
		SORT_ORDER=coordinate

	# Remove duplicate reads generated during sequencing
	java -jar /Path/To/picard/picard.jar MarkDuplicates \
		INPUT=${SAMPLE}_Trimmed_P_bwa_aligned_sorted.bam \
		OUTPUT=${SAMPLE}_Trimmed_aligned_sorted_merged_dedup.bam \
		METRICS_FILE=${SAMPLE}_dedup_metrics.txt

	# Indexing the final BAM file
	java -jar /Path/To/picard/picard.jar BuildBamIndex \
		INPUT=${SAMPLE}_Trimmed_aligned_sorted_merged_dedup.bam

done
