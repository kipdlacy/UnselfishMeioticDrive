#! /usr/bin/bash

# This file contains the commands I used to call variants among short read whole genomes that had previously been aligned and processed.

GATK=/home/kip/tools/gatk-4.2.0.0/gatk

# Call variants using HaplotypeCaller from GATK 4.2
# Note that this process goes much faster if you paralelize by chromosome
# For simplicity, I don't include the names of all samples used in this study, and instead just mention hypothetical samples.
# $GATK HaplotypeCaller \
# 	-R /Path/To/Reference/Genome/Obir.assembly.v5.4.fasta \
# 	-I Sample1_Trimmed_aligned_sorted_merged_dedup.bam \
# 	-I Sample2_Trimmed_aligned_sorted_merged_dedup.bam \
# 	-I Sample3_Trimmed_aligned_sorted_merged_dedup.bam \
# 	-O AllShortReadGenomes.vcf

# These results also need to be filtered before analysis.
# The specific filtering steps I performed for each analysis can be found the appropriate directory.
# But I always performed the following basic filtering steps.

# Select only single nucleotide polymorphisms, because these are the most reliably identified variant with short-read sequencing data
$GATK SelectVariants \
	-V ShortReadGenomes.vcf \
	--select-type-to-include SNP \
	-O ShortReadGenomes_SNPs.vcf

# Filter the variants according to the GATK best practice recommendations
# See this link for details: "https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
$GATK VariantFiltration \
	-V ShortReadGenomes_SNPs.vcf \
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
	-O ShortReadGenomes_SNPs_GATKfilt.vcf

# Remove filtered sites from the file so that only unfiltered sites are used for down stream processing.
$GATK SelectVariants \
	-V ShortReadGenomes_SNPs_GATKfilt.vcf \
	--exclude-filtered \
	-O ShortReadGenomes_SNPs_GATKfilt_RmFilt.vcf

# The O. biroi reference genome is overall quite good, but contains assembly errors in some places. Variants found in these regions
# cannot be trusted, and need to be excluded from analyses. Of particular concern for this study are falsely heterozygous sites, which
# occur in regions of the reference genome that were formed by incorrect assembly of two distinct genomic regions with variation as
# a single region. The variants between these two regions will appear to be heterozygous in all samples. In our study, we were able
# to rule out many of these false variants by using haploid male genomes as an internal control. Haploids should have no heterozygosity
# in their nuclear genomes. By identifying sites that were found to be heterozygous in haploid samples, we could eliminate those
# sites from downstream analyses.

# This process was performed by using GATK SelectVariants to include only haploid male samples.
$GATK SelectVariants \
	-V ShortReadGenomes_SNPs_GATKfilt_RmFilt.vcf \
	--sample-name SM03 \
	--sample-name SM06 \
	--sample-name SM65 \
	--sample-name SM74 \
	--sample-name SM75 \
	--sample-name SM76 \
	--sample-name SM83 \
	--sample-name SM87 \
	-O ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf

# Then, we write these genotypes to a VCF-style table
$GATK VariantsToTable \
	-V ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf \
	-F CHROM -F POS -F REF -F ALT -F FILTER \
	-GF GT \
	-O ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf.table

# Then we identify all sites at which at least one haploid male is heterozygous, and write those to a new file.
python ./KeepOnlySitesWhereAtLeastOneSampleIsHeterozygous.py ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf.table

# Then we convert this into a simple table that just shows the chromosomal position of each of these sites.
python ./Process_VCFtoPosScreener.py ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf.table_AtLeastOneSampleHet

# We can then use vcftools to screen out these sites from the VCF we want to perform downstream analyses on.
vcftools \
	--vcf ShortReadGenomes_SNPs_GATKfilt_RmFilt.vcf \
	--exclude-positions ShortReadGenomes_SNPs_GATKfilt_RmFilt_HapMalesOnly.vcf.table_AtLeastOneSampleHet.PosScreener \
	--recode \
	--recode-INFO-all \
	--out ShortReadGenomes_SNPs_GATKfilt_RmFilt_NoHapHet

# From here, the file "ShortReadGenomes_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf" excludes all sites heterozygous in haploid males.

# Removing all files except for this final file to clear up space.
#find ./ -name "ShortReadGenomes_*" ! -name "*NoHapHet.recode.vcf" -exec rm {} \;
