#! /usr/bin/bash

# This file contains the commands I used to call variants among short read whole genomes that had previously been aligned and processed.
# For simplicity, I don't include the names of all samples used in this study, and instead loop across hypothetical samples.

# Call variants using HaplotypeCaller from GATK 4.2
# Note that this process goes much faster if you paralelize by chromosome
/Path/To/gatk-4.2.0.0/gatk HaplotypeCaller \
	-R /Path/To/Reference/Genome/Obir.assembly.v5.4.fasta \
	-I Sample1_Trimmed_aligned_sorted_merged_dedup.bam \
	-I Sample2_Trimmed_aligned_sorted_merged_dedup.bam \
	-I Sample3_Trimmed_aligned_sorted_merged_dedup.bam \
	-O AllSamples.vcf

# These results also need to be filtered before analysis.
# The specific filtering steps I performed for each analysis can be found the appropriate directory.
# But I always performed the following basic filtering steps.

# Select only single nucleotide polymorphisms, because these are the most reliably identified variant with short-read sequencing data
/Path/To/gatk-4.2.0.0/gatk SelectVariants \
	-V AllSamples.vcf \
	--select-type-to-include SNP \
	-O AllSamples_SNPs.vcf

# Filter the variants according to the GATK best practice recommendations
# See this link for details: "https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
/Path/To/gatk-4.2.0.0/gatk VariantFiltration \
	-V AllSamples_SNPs.vcf \
	--filter-name "QD" \
	--filter-expression "QD < 2.0" \
	--filter-name "FS" \
	--filter-expression "FS > 55.0" \
	--filter-name "MQ" \
	--filter-expression "MQ < 55.0" \
	--filter-name "MQRankSum" \
	--filter-expression "MQRankSum < -8.5" \
	--filter-name "ReadPosRankSum" \
	--filter-expression "ReadPosRankSum < -6.0" \
	-O AllSamples_SNPs_GATKfilters.vcf

# Remove filtered sites from the file so that only unfiltered sites are used for down stream processing. 
/Path/To/gatk-4.2.0.0/gatk SelectVariants \
	-V AllSamples_SNPs_GATKfilters.vcf \
	--exclude-filtered \
	-O AllSamples_SNPs_GATKfilters_RmFilt.vcf
