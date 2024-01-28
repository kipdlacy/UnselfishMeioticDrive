#! /usr/bin/bash

### 20220413
### This file is designed to remove known problematic SNPs from phased VCFs

GATK=/home/kip/tools/gatk-4.2.0.0/gatk

### need to start by copying all the phased.VCF files to the appropriate directories.
### Grab data and make folders
### we are using internal codes
Samples=("FEB_T505" "FEB_T506" "FEB_T507" "FEB_T508" "FEB_T509" "FEB_T510" "FEB_T511" "FEB_T512" "FEB_T513" "FEB_T514" "FEB_T515" "FEB_T516")
# Loop through each SAMPLE
for SAMPLE in "${Samples[@]}"; do
	# Extract the part after "FEB_" from the SAMPLE name
	SAMPLESuffix=${SAMPLE#FEB_}
	# Create a directory for the SAMPLE in the ./data/Samples/ directory
	mkdir -p "./data/Samples/${SAMPLE}/"
	# Copy VCF files from the source directory to the newly created SAMPLE directory
	cp "/path/to/output/directoryFEB/results/tellsortOUT/${SAMPLESuffix}/${SAMPLESuffix}_temp/step2_run_hapcut2_ustpipeline/s3_hapcut_output/hapblock_chr*.phased.VCF" "./data/Samples/${SAMPLE}/"
done
Samples=("NOV_T501" "NOV_T503" "NOV_T504" "NOV_T505" "NOV_T506")
# Loop through each SAMPLE
for SAMPLE in "${Samples[@]}"; do
	# Extract the part after "FEB_" from the SAMPLE name
	SAMPLESuffix=${SAMPLE#NOV_}
	# Create a directory for the SAMPLE in the ./data/Samples/ directory
	mkdir -p "./data/Samples/${SAMPLE}/"
	# Copy VCF files from the source directory to the newly created SAMPLE directory
	cp "/path/to/output/directoryNOV/results/tellsortOUT/${SAMPLESuffix}/${SAMPLESuffix}_temp/step2_run_hapcut2_ustpipeline/s3_hapcut_output/hapblock_chr*.phased.VCF" "./data/Samples/${SAMPLE}/"
done

### For each linked read genome (here showing internal codes)
for i in NOV_T501 NOV_T503 NOV_T504 NOV_T505 NOV_T506 FEB_T505 FEB_T506 FEB_T507 FEB_T508 FEB_T509 FEB_T510 FEB_T511 FEB_T512 FEB_T513 FEB_T514 FEB_T515 FEB_T516; do
	### Separately for each chromosome
	for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
		### Remove all sites that overlap with predicted transposable elements (TEs)
		### Such sites are often erroneously assembled in reference genomes, and therefore variants that fall within them
		### are likely to represent sequences from different parts of the reference genome rather than allelic variation
		### therefore, we cannot use such sites to identify crossover recombination events by aligning short reads to those sites
		vcftools --vcf ./data/Samples/$i/hapblock_chr$CHROM.phased.VCF \
			--exclude-bed ../DataFiles/data/Obir.assembly.v5.4.fasta.out_TE.bed \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$CHROM.phased_TEfiltered
		### Then identify all those SNPs that are present within predicted TEs
		vcftools --vcf ./data/Samples/$i/hapblock_chr$CHROM.phased.VCF \
			--bed ../DataFiles/data/Obir.assembly.v5.4.fasta.out_TE.bed \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$CHROM.phased_OnlyTE
		### Then remove all sites at which any haploid has been determined to be heterozygous
		vcftools --vcf ./data/Samples/$i/hapblock_chr$CHROM.phased_TEfiltered.recode.vcf \
			--exclude-positions ../DataFiles/data/HaploidMaleHet.PosScreener \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$CHROM.phased_TEfiltered.recode_NoHapHet
		### Then identify all such positions
		vcftools --vcf ./data/Samples/$i/hapblock_chr$CHROM.phased_TEfiltered.recode.vcf \
			--positions ../DataFiles/data/HaploidMaleHet.PosScreener \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$CHROM.phased_TEfiltered.recode_OnlyHapHet
	done
done

### In a new loop across chromosomes for each sample filter by depth.
### We have to each sample separately, because each sample has a different genome-wide read depth avg.
### We always want to filter out variants with read depths below 20, as these are likely poor quality variants that will not be trustworthy for
### heterozygosity calling. We want to filter out sites that are above 2x the genome-wide average for each sample
declare -A maxDP=([NOV_T501]=141.9158 [NOV_T503]=348.096 [NOV_T504]=385.406 [NOV_T505]=165.7718 [NOV_T506]=130.7658 [FEB_T505]=142.9356 [FEB_T506]=97.5684 [FEB_T507]=192.2282 [FEB_T508]=173.817 [FEB_T509]=172.5012 [FEB_T510]=155.458 [FEB_T511]=145.1288 [FEB_T512]=135.8142 [FEB_T513]=149.2908 [FEB_T514]=135.1938 [FEB_T515]=173.3138 [FEB_T516]=158.2244)
### Loop over chromosomes
for CHROM in {1..14}; do
	### Loop over each sample in the array
	for SAMPLE in "${!maxDP[@]}"; do
		### Define file paths
		INPUT="./data/Samples/${SAMPLE}/hapblock_chr${CHROM}.phased_TEfiltered.recode_NoHapHet.recode.vcf"
		OUTPUT="./data/Samples/${SAMPLE}/hapblock_chr${CHROM}.phased_TEfiltered.recode_NoHapHet.recode"
		### Removing all sites below 20 DP and above 2x Mean Depth
		vcftools --vcf "$INPUT" --min-meanDP 20.0 --max-meanDP "${maxDP[$SAMPLE]}" --recode --recode-INFO-all --out "${OUTPUT}_NonAbove2xMeanDepth"
		### Then producing files with just those values
		vcftools --vcf "$INPUT" --max-meanDP 19.99999 --recode --recode-INFO-all --out "${OUTPUT}_Max20Depth"
		vcftools --vcf "$INPUT" --min-meanDP "${maxDP[$SAMPLE]}" --recode --recode-INFO-all --out "${OUTPUT}_Min2xMeanDepth"
	done
done

### FOR EACH SAMPLE, CREATE A FILE WITH ONLY THE BAD SITES
### This will be used for screening putative crossover events in downstream analyses
for i in NOV_T501 NOV_T503 NOV_T504 NOV_T505 NOV_T506 FEB_T505 FEB_T506 FEB_T507 FEB_T508 FEB_T509 FEB_T510 FEB_T511 FEB_T512 FEB_T513 FEB_T514 FEB_T515 FEB_T516; do
	### Looping through chromosomes
	for cHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
		### Create just a simple file that only contains the name of the chromosome and the position number of the variant
		echo "CHROM\tPOS" > ./data/Samples/$i/chr$cHROM'_AllScreenedSitesUnordered_NoIndels.txt'
		### Add sites present in TEs, add sites found to be heterozygous in haploids, add sites with read depth < 20, add sites with read depth > the genome-wide mean
		for FILE in _OnlyTE.recode.vcf _TEfiltered.recode_OnlyHapHet.recode.vcf _TEfiltered.recode_NoHapHet.recode_Max20Depth.recode.vcf _TEfiltered.recode_NoHapHet.recode_Min2xMeanDepth.recode.vcf; do
			### Grab only SNPs
			$GATK SelectVariants \
				-V ./data/Samples/$i/hapblock_chr$cHROM.phased$FILE \
				--select-type-to-include SNP \
				-O ./data/Samples/$i/hapblock_chr$cHROM.phased$FILE'_SNPsOnly.vcf'
			### Write a table
			$GATK VariantsToTable \
				-V ./data/Samples/$i/hapblock_chr$cHROM.phased$FILE'_SNPsOnly.vcf' \
				-F CHROM -F POS \
				-O ./data/Samples/$i/hapblock_chr$cHROM.phased$FILE'_SNPsOnly.vcf'.table
			### then add all sites from that table.
			tail -n +2 ./data/Samples/$i/hapblock_chr$cHROM.phased$FILE'_SNPsOnly.vcf'.table >> ./data/Samples/$i/chr$cHROM'_AllScreenedSitesUnordered_NoIndels.txt'
		done
		#### remove indels from focal file, and create a table version with GQ and PQ
		### (Genotype quality and phasing quality)
		$GATK SelectVariants \
			-V ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf \
			--select-type-to-include SNP \
			-O ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'
		$GATK VariantsToTable \
			-V ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf' \
			-F CHROM -F POS \
			-GF GQ -GF PQ \
			-O ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ
		### Custom script to filter based on these values
		### I identify everything that's less than perfect
		python ./src/FilterByFieldFromVCFtable.py ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ GQ LessThan 99
		python ./src/FilterByFieldFromVCFtable.py ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ PQ LessThan 100
		### Then filter it out using vcftools
		vcftools --vcf ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf' \
			--exclude-positions ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ_GQLessThan99 \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99
		vcftools --vcf ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode.vcf \
			--exclude-positions ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ_PQLessThan100 \
			--recode --recode-INFO-all \
			--out ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100
		### Write those Sites to a table
		$GATK VariantsToTable \
			-V ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
			-F CHROM -F POS \
			-O ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table
		### And then put them in the ~allbadvariants~ file
		tail -n +2 ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ_GQLessThan99 >> ./data/Samples/$i/chr$cHROM'_AllScreenedSitesUnordered_NoIndels.txt'
		tail -n +2 ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode.vcf'_SNPsOnly.vcf'.table_GQPQ_PQLessThan100 >> ./data/Samples/$i/chr$cHROM'_AllScreenedSitesUnordered_NoIndels.txt'
		### Separately, write allelic depths to a table for a different downstream screening tool
		$GATK VariantsToTable \
			-V ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
			-F CHROM -F POS \
			-GF AD \
			-O ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table_AD
		### Turn that table into a table of Minor Allelic frequencies (MAF) using a custom script
		python ./src/TurnADTableToTableOfMAF.py ./data/Samples/$i/hapblock_chr$cHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table_AD
	done
done

### Generate lists of sites that have zero read depth, or read depth greater than the genome-wide mean for each sample
### These will be used in downstream screening of putative recombination events
### We have separate loops for the "FEB" and "NOV" series samples, because the file structure output from the UST pipelines
### produces folders with the name of the barcode only, and thus, overlap in names between the two series

# Define an associative array with the maximum depth values for each sample
declare -A NOVmaxDP=([NOV_T501]=141.9158 [NOV_T503]=348.096 [NOV_T504]=385.406 [NOV_T505]=165.7718 [NOV_T506]=130.7658)
# Loop through each sample in the maxDepths associative array
# SAMPLE variable holds the current SAMPLE name (e.g., FEB_T505)
# maxDepths[$SAMPLE] accesses the maximum depth value for the current SAMPLE
for SAMPLE in "${!NOVmaxDP[@]}"; do
	# Constructing the path to the BAM file for the current SAMPLE
	SAMPLEsuffix=${SAMPLE#NOV_}
	BamFilePath="/path/to/output/directoryNOV/results/tellsortOUT/${SAMPLEsuffix}/${SAMPLEsuffix}_temp/phased_sorted.bam"
	# File paths for storing all depth data and quality-filtered depth data
	DepthFile="./data/Samples/${SAMPLE}/AllDepth.txt"
	QualDepthFile="./data/Samples/${SAMPLE}/AllQualDepth.txt"
	# Run samtools to get depth information for all positions in the BAM file
	nohup samtools depth -a "${BamFilePath}" >> "${DepthFile}"
	# Run Python script to extract zero depth sites
	python ./src/ExtractZeroDepthSitesFromAllDepth.py "${DepthFile}"
	# Run Python script to extract sites with read depth greater than 2x the genome-wide mean
	python ./src/ExtractHighDepthSitesFromAllDepth.py "${DepthFile}" "${maxDepths[$SAMPLE]}"
	# Remove the depth file as it is no longer needed
	rm "${DepthFile}"
	# Repeat the process for quality-filtered depth information
	nohup samtools depth -a -Q 1 "${BamFilePath}" >> "${QualDepthFile}"
	python ./src/ExtractZeroDepthSitesFromAllDepth.py "${QualDepthFile}"
	rm "${QualDepthFile}"
done

# Define an associative array with the maximum depth values for each sample
declare -A FEBmaxDP=([FEB_T505]=142.9356 [FEB_T506]=97.5684 [FEB_T507]=192.2282 [FEB_T508]=173.817 [FEB_T509]=172.5012 [FEB_T510]=155.458 [FEB_T511]=145.1288 [FEB_T512]=135.8142 [FEB_T513]=149.2908 [FEB_T514]=135.1938 [FEB_T515]=173.3138 [FEB_T516]=158.2244)
# Loop through each sample in the maxDepths associative array
# SAMPLE variable holds the current SAMPLE name (e.g., FEB_T505)
# maxDepths[$SAMPLE] accesses the maximum depth value for the current SAMPLE
for SAMPLE in "${!FEBmaxDP[@]}"; do
	# Constructing the path to the BAM file for the current SAMPLE
	SAMPLEsuffix=${SAMPLE#FEB_}
	BamFilePath="/path/to/output/directoryFEB/results/tellsortOUT/${SAMPLEsuffix}/${SAMPLEsuffix}_temp/phased_sorted.bam"
	# File paths for storing all depth data and quality-filtered depth data
	DepthFile="./data/Samples/${SAMPLE}/AllDepth.txt"
	QualDepthFile="./data/Samples/${SAMPLE}/AllQualDepth.txt"
	# Run samtools to get depth information for all positions in the BAM file
	nohup samtools depth -a "${BamFilePath}" >> "${DepthFile}"
	# Run Python script to extract zero depth sites
	python ./src/ExtractZeroDepthSitesFromAllDepth.py "${DepthFile}"
	# Run Python script to extract sites with read depth greater than 2x the genome-wide mean
	python ./src/ExtractHighDepthSitesFromAllDepth.py "${DepthFile}" "${maxDepths[$SAMPLE]}"
	# Remove the depth file as it is no longer needed
	rm "${DepthFile}"
	# Repeat the process for quality-filtered depth information
	nohup samtools depth -a -Q 1 "${BamFilePath}" >> "${QualDepthFile}"
	python ./src/ExtractZeroDepthSitesFromAllDepth.py "${QualDepthFile}"
	rm "${QualDepthFile}"
done
