#! /usr/bin/bash

# This script is designed to perform pairwise analysis of changes in heterozygosity between genomes,
# as well as to produce files that report the proportion of 10kb windows that have at least one heterozygous SNP

# Telling bash where to find the GATK software
GATK=/home/kip/tools/gatk-4.2.0.0/gatk

# # Creating folders to store data files and results from the analysis
mkdir ./data/
cp ../DataFiles/data/10000bpWindowsNonOverlapping.txt ./data/

# Select samples for this analysis from a file that contains all variants used in study
$GATK SelectVariants \
	-V ../DataFiles/data/ShortReadGenomes.vcf \
	--sample-name MD-B-1-D \
	--sample-name MD-B-1-M \
	--sample-name MD-B-7-D1 \
	--sample-name MD-B-7-M \
	-O ./data/ShortRead_MDpairs.vcf

# Select only single nucleotide polymorphisms, because these are the most reliably identified variant with short-read sequencing data
$GATK SelectVariants \
	-V ./data/ShortRead_MDpairs.vcf \
	--select-type-to-include SNP \
	-O ./data/ShortRead_MDpairs_SNPs.vcf

# Filter the variants according to the GATK best practice recommendations
# See this link for details: "https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
$GATK VariantFiltration \
	-V ./data/ShortRead_MDpairs_SNPs.vcf \
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
	-O ./data/ShortRead_MDpairs_SNPs_GATKfilt.vcf

# Remove filtered sites from the file so that only unfiltered sites are used for down stream processing.
$GATK SelectVariants \
	-V ./data/ShortRead_MDpairs_SNPs_GATKfilt.vcf \
	--exclude-filtered \
	-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt.vcf

# The O. biroi reference genome is overall quite good, but contains assembly errors in some places. Variants found in these regions
# cannot be trusted, and need to be excluded from analyses. Of particular concern for this study are falsely heterozygous sites, which
# occur in regions of the reference genome that were formed by incorrect assembly of two distinct genomic regions with variation as
# a single region. The variants between these two regions will appear to be heterozygous in all samples. In our study, we were able
# to rule out many of these false variants by using haploid male genomes as an internal control. Haploids should have no heterozygosity
# in their nuclear genomes. By identifying sites that were found to be heterozygous in haploid samples, we could eliminate those
# sites from downstream analyses.

# We can then use vcftools to screen out these sites from the VCF we want to perform downstream analyses on.
# From here, the file "ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf" excludes all sites heterozygous in haploid males.
vcftools \
	--vcf ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt.vcf \
	--exclude-positions ../DataFiles/data/HaploidMaleHet.PosScreener \
	--recode \
	--recode-INFO-all \
	--out ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet

# Removing all files except for this final file to clear up space.
find ./ -name "./data/ShortRead_MDpairs_*" ! -name "*NoHapHet.recode.vcf" -exec rm {} \;

# # For analyses of heterozygosity in parthenogenetic lineages, it is important to only consider sites that were putatively ancestrally heterozygous.
# # In this case, we mean variants for which at least one sample is heterozygous, and all other samples are either also heterozygous (with an identical genotype)
# # or are homozygous for one of the alleles included in the genotype of the heterozygous sample(s). Such variants would have been heterozygous in the common ancestor of all samples
# # and may have lost heterozygosity in one or more samples. Note that this method will also identify some sites that have gained heterozygosity via point mutation since
# # the foundation of the parthenogenetic line. We then must also assure rigor in these genotype calls by screening out putatively heterozygous or homozygous calls of poor quality.

# First we write a genotype table ("GT" field) for use in identifying "putatively ancestrally heterozygous" variants
$GATK VariantsToTable \
	-V ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf \
	-F CHROM -F POS -F REF -F ALT -F FILTER \
	-GF GT \
	-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table
# 
# # Then we write an Allelic Depth table ("AD" field) for use in ruling out sites with poor quality heterozygosity or homozygosity calls
$GATK VariantsToTable \
	-V ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf \
	-F CHROM -F POS -F REF -F ALT -F FILTER \
	-GF AD \
	-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table
# 
# # Using a custom python script we produce a new output file that contains only variants that are putatively ancestrally heterozygous
python CleanData_KeepPutativelyAncestralHet.py ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table
# 
# # Screen sites based on allelic depth to remove variants with questionable heterozygosity or homozygosity calls
python ScreenSNPsonADnDP.py ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table

# Convert the resulting table of high quality SNPs to a bed style "positions" file to screen against
python ../AlignmentAndVariantCallingWGS/Process_VCFtoPosScreener.py  ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs

# Use vcftools to create a new VCF file that contains only those high quality, trusted positions
vcftools \
	--vcf ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf \
	--positions ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener \
	--recode \
	--recode-INFO-all \
	--out ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs

# Now that we have files with trustworthy variants that are putatively ancestrally heterozygous, we can analyze heterozygosity patterns across individuals.
# We will start by writing separate VCF files for each pairwise comparison, and writing the genotype data to a table. 
# Then we will find all differences in heterozygosity between the individuals in each pair.

# Initialize an array to hold unique sample names
declare -A UniqueSamples
# List of comparisons
Comparisons="MD-B-1-M_v_MD-B-1-D MD-B-7-M_v_MD-B-7-D1 MD-B-1-M_v_MD-B-7-M"
# Loop through each comparison
for COMP in $Comparisons; do
	# Split the comparison into two samples using '_v_' as delimiter
	Sample1=${COMP%_v_*}
	Sample2=${COMP#*_v_}
	# Store the samples in the associative array to ensure uniqueness
	UniqueSamples["$Sample1"]=1
	UniqueSamples["$Sample2"]=1
	$GATK SelectVariants \
		-V ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf \
		--sample-name $Sample1 \
		--sample-name $Sample2 \
		--exclude-non-variants \
		--remove-unused-alternates \
		-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_$COMP'.vcf'
	#Write variants to table
	$GATK VariantsToTable \
		-V ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_$COMP'.vcf' \
		-F CHROM -F POS -F REF -F ALT -F FILTER \
		-GF GT \
		-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_$COMP'.vcf.table'
	#Note that the output of this file "_Diffs" will still contain false positives that need to be manually screened by visual inspection in IGCV
	python Process_GetDiffs.py ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_$COMP'.vcf''.table'
done

# Okay, now we want to summarize the proportion of the genome that contains heterozygosity for each sample
# To do this, we will first obtain a table of high quality variants for all samples, 
# Then, we will write files for each individual that contain all sites at which that sample is heterozygous,
# And we will count the number of adjacent (i.e., non overlapping) windows that contain at least one heterozygous SNP.

# Writing all samples to a table
$GATK VariantsToTable \
	-V ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf \
	-F CHROM -F POS -F REF -F ALT -F FILTER \
	-GF GT \
	-O ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table

# Writing one file for each individual that contains all sites at which each individual is heterozygous
python WriteIndividualSampleHetFiles.py ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table

# For each sample, tallying the number of heterozygous sites and recording a summary
echo -e "SampleName\tWindowSize(bp)\tWindowsHet\tTotalWindows\tProportionWindowsHet" > ./data/Summary_ShortRead_10000bpWindowsHet.txt
for SAMPLE in "${!UniqueSamples[@]}"; do
	python Analysis_PropGenomeHet.py ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_$SAMPLE'Het' ./data/10000bpWindowsNonOverlapping.txt
	tail -n +2 ./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_$SAMPLE'Het_Prop_10000bpWindowsHet' >> ./data/Summary_ShortRead_10000bpWindowsHet.txt
	echo >> ./data/Summary_ShortRead_10000bpWindowsHet.txt
done
