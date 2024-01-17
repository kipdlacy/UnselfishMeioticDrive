#! /usr/bin/bash

### This script is designed to clean and filter VCF files for analyses of haplotypic recombination among individuals of unknown pedigree.
### For haploid males, haplotypic recombination can be directly inferred, while for diploid females, haplotypic recombination can be inferred from LOH

GATK=/home/kip/tools/gatk-4.2.0.0/gatk

# mkdir ./data/

# #######################################################################
# ###
# ###
# ###  First we will perform within-colony analyses of recombination
# ###
# ###
# #######################################################################
# 
# ### We start by obtaining working files that contain just the samples we need for each analysis. Note that these use internal sample codes
# ### One analysis will be diploid females and haploid males from Line A colony C16
# $GATK SelectVariants -V ../DataFiles/data/ShortReadGenomes.vcf --sample-name C16B-3-1 --sample-name C16B-4-2 --sample-name SM65 --sample-name SM76 --sample-name SM83 --sample-name SM87 --sample-name W06 --sample-name W07 --select-type-to-include SNP -O ./data/C16_SNPs.vcf
# ### Another will be diploid females and haploid males from Line B colony STC6
# $GATK SelectVariants -V ../DataFiles/data/ShortReadGenomes.vcf --sample-name SM03 --sample-name SM74 --sample-name SM75 --sample-name STC6L --sample-name W04 --sample-name W05 --select-type-to-include SNP -O ./data/STC6_SNPs.vcf
# 
# for COL in C16 STC6; do
# 	# Filter according to GATK best practices
# 	$GATK VariantFiltration -V ./data/$COL'_SNPs.vcf' --filter-name "QD" --filter-expression "QD < 2.0" --filter-name "FS" --filter-expression "FS > 55.0" --filter-name "MQ" --filter-expression "MQ < 55.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -8.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -6.0" -O ./data/$COL'_SNPs_GATKfilters.vcf'
# 	# Remove filtered sites
# 	$GATK SelectVariants -V ./data/$COL'_SNPs_GATKfilters.vcf' --exclude-filtered -O ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf'
# 	# Exclude sites where haploid samples are heterozygous (erroneously assembled regions of the reference genome)
# 	vcftools --vcf ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf' --exclude-positions ../DataFiles/data/HaploidMaleHet.PosScreener  --recode --recode-INFO-all --out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet'
# 	# Write genotypes to a table
# 	$GATK VariantsToTable -V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF GT -O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# 	# Write allelic depths to a table
# 	$GATK VariantsToTable -V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF AD -O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
# 	# Keep only sites that are putatively ancestrally heterozygous for the parthenogenetic lineage	
# 	python ./CleanData_KeepPutativelyAncestralHet_ModForHaploids.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# 	# Remove null variant calls (./. or ./C), and filter based on allelic depth
# 	python ../LossOfHeterozygosity_KnownPedigreePairs/ScreenSNPsonADnDP.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet' ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
# 	# Write a simple file to screen the VCF by
# 	python ../AlignmentAndVariantCallingWGS/Process_VCFtoPosScreener.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs'
# 	# Screen the VCF
# 	vcftools --vcf ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' --positions ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener' --recode --recode-INFO-all --out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs'
# done
# 
# ## Now that the files are filtered, write files with relevant samples for analysis
# # Create a file with all Line A C16 Diploid Females
# $GATK SelectVariants -V ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name C16B-3-1 --sample-name C16B-4-2 --sample-name W06 --sample-name W07 --exclude-non-variants --remove-unused-alternates -O ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf
# # Create a file with all Line A C16 Haploid Males
# $GATK SelectVariants -V ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name SM65 --sample-name SM76 --sample-name SM83 --sample-name SM87 --exclude-non-variants --remove-unused-alternates -O ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf
# # Create a file with all Line B STC6 Diploid Females
# $GATK SelectVariants -V ./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name STC6L --sample-name W04 --sample-name W05 -O ./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf
# # Create a file with all Line B STC6 Haploid Males
# $GATK SelectVariants -V ./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name SM03 --sample-name SM74 --sample-name SM75 -O ./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf
# 
# for COL in C16 STC6; do
# 	## Prepare data for analysis
# 	for TYPE in DipFem HapMale; do
# 		# write genotype tables
# 		$GATK VariantsToTable -V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF GT -O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table'
# 		# Add contig names to the tables. The contig names help to determine when the start and endpoints of haplotype blocks or runs of homozygosity occur.
# 		python ./CleanData_AddContigNames.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table' ../UsefulFiles/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# 	done
# 	## find runs of homozygosity among the group of diploids
# 	python ./LOH_GroupOfDiploidsHomRunFinder.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs' ./data/
# 	## create a file for HaplotypeBlock presentation
# 	python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt'
# 	## find crossovers among haploid male genomes
# 	python ./Recomb_HapPhaseSwitchScanner.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs'
# 	## using a realistic cutoff for Noncrossover (NCO) tract length (1500) and a conservative one (5000)
# 	for NCOcutoff in 1500 5000; do
# 		###### First for the diploids
# 		## create a file that removes all NCO-associated LOH and keeps only putative crossover (CO)-associated LOH
# 		python ./PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' $NCOcutoff
# 		## create a file that removes all CO-associated LOH and keeps only putative NCO-associated LOH
# 		python ./PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' $NCOcutoff
# 		###### then for the haploids
# 		## create a file that removes all NCO-associated LOH and keeps only putative crossover (CO)-associated LOH
# 		python ./PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' $NCOcutoff
# 		## create a file that removes all CO-associated LOH and keeps only putative NCO-associated LOH
# 		python ./PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' $NCOcutoff
# 		## Once you've removed NCOs, there will be some CO-associated haplotype blocks that were previously interrupted by NCOs. These need to be merged.
# 		python ./PhaseClean_MergeAdjacents.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP.txt'
# 		## Once you've done that merging, we can tally the number of recombination events detected.
# 		python ./PhaseClean_CountPhaseSwitches.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents.txt'
# 		## create a file for HaplotypeBlock presentation
# 		python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan'$NCOcutoff'BP.txt'
# 		python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents.txt'
# 		## create a file that contains all recombination events, including crossovers and gene conversions (NCOs)
# 		cat ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt' > ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# 		echo '' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# 		tail -n +2 ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan'$NCOcutoff'BP_FirstBlocksOnContigRemoved.txt' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# 	done
# done
# ### Note that these files need to be manually curated to remove false positives, such as long LOH tracts that are only supported by a few SNPs.


######################################################################################################
###
###
###  Now we will perform analyses of recombination across multiple colonies for each clonal line
###
###
######################################################################################################


$GATK SelectVariants -V /home/antqueen/genomics/experiments/analyses/KDL20211215_VarCall4AllIndieGenomesToDateGATK4.2/results/AllIndidualsAsOf_20211215_AllChroms_SNPs.vcf --sample-name C16B-3-1 --sample-name C16B-4-2 --sample-name W06 --sample-name W07 --sample-name C18 --sample-name C1A-1 --sample-name SM65 --sample-name SM76 --sample-name SM83 --sample-name SM87 --sample-name KL46 --sample-name KL66 -O ./data/LineAmultiCol_SNPs.vcf

$GATK SelectVariants -V /home/antqueen/genomics/experiments/analyses/KDL20211215_VarCall4AllIndieGenomesToDateGATK4.2/results/AllIndidualsAsOf_20211215_AllChroms_SNPs.vcf --sample-name W04 --sample-name W05 --sample-name STC6L --sample-name STC1_W --sample-name STC12_W --sample-name SM03 --sample-name SM74 --sample-name SM75 --sample-name KL70 --sample-name SM56 -O ./data/LineBmultiCol_SNPs.vcf

for LINE in LineAmultiCol LineBmultiCol; do
	# Filter according to GATK best practices
	$GATK VariantFiltration -V ./data/$LINE'_SNPs.vcf' --filter-name "QD" --filter-expression "QD < 2.0" --filter-name "FS" --filter-expression "FS > 55.0" --filter-name "MQ" --filter-expression "MQ < 55.0" --filter-name "MQRankSum" --filter-expression "MQRankSum < -8.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -6.0" -O ./data/$LINE'_SNPs_GATKfilters.vcf'
	# Remove filtered sites
	$GATK SelectVariants -V ./data/$LINE'_SNPs_GATKfilters.vcf' --exclude-filtered -O ./data/$LINE'_SNPs_GATKfilt_RmFilt.vcf'
	# Exclude sites where haploid samples are heterozygous (erroneously assembled regions of the reference genome)
	vcftools --vcf ./data/$LINE'_SNPs_GATKfilt_RmFilt.vcf' --exclude-positions ../DataFiles/data/HaploidMaleHet.PosScreener  --recode --recode-INFO-all --out ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet'
	# Write genotypes to a table
	$GATK VariantsToTable -V ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF GT -O ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
	# Write allelic depths to a table
	$GATK VariantsToTable -V ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF AD -O ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
	# Keep only sites that are putatively ancestrally heterozygous for the parthenogenetic lineage	
	python ./CleanData_KeepPutativelyAncestralHet_ModForHaploids.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
	# Remove null variant calls (./. or ./C), and filter based on allelic depth
	python ../LossOfHeterozygosity_KnownPedigreePairs/ScreenSNPsonADnDP.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet' ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
	# Write a simple file to screen the VCF by
	python ../AlignmentAndVariantCallingWGS/Process_VCFtoPosScreener.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs'
	# Screen the VCF
	vcftools --vcf ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' --positions ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener' --recode --recode-INFO-all --out ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs'
done

# ## Now that the files are filtered, write files with relevant samples for analysis
# # Create a file with all Line A Diploid Females
$GATK SelectVariants -V ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name C16B-3-1 --sample-name C16B-4-2 --sample-name W06 --sample-name W07 --sample-name C18 --sample-name C1A-1 --exclude-non-variants --remove-unused-alternates -O ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf
# # Create a file with all Line B Diploid Females
$GATK SelectVariants -V ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name STC6L --sample-name W04 --sample-name W05 --sample-name STC1_W --sample-name STC12_W --exclude-non-variants --remove-unused-alternates -O ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf
# # Create a file with all Line A Haploid Males
$GATK SelectVariants -V ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name SM65 --sample-name SM76 --sample-name SM83 --sample-name SM87 --sample-name KL46 --sample-name KL66 --exclude-non-variants --remove-unused-alternates -O ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf
# # Create a file with all Line B Haploid Males
$GATK SelectVariants -V ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf --sample-name SM03 --sample-name SM74 --sample-name SM75 --sample-name KL70 --sample-name SM56 --exclude-non-variants --remove-unused-alternates -O ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf


for LINE in LineAmultiCol LineBmultiCol; do
	## Prepare data for analysis
	for TYPE in DipFem HapMale; do
		# write genotype tables
		$GATK VariantsToTable -V ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf' -F CHROM -F POS -F REF -F ALT -F FILTER -GF GT -O ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table'
		# Add contig names to the tables. The contig names help to determine when the start and endpoints of haplotype blocks or runs of homozygosity occur.
		python ./CleanData_AddContigNames.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table' ../UsefulFiles/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
	done
	## find runs of homozygosity among the group of diploids
	python ./LOH_GroupOfDiploidsHomRunFinder.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs' ./data/
	## create a file for HaplotypeBlock presentation
	python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt'
	## find crossovers among haploid male genomes
	python ./Recomb_HapPhaseSwitchScanner.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs'
	## using a realistic cutoff for Noncrossover (NCO) tract length (1500) and a conservative one (5000)
	for NCOcutoff in 1500 5000; do
		###### First for the diploids
		## create a file that removes all NCO-associated LOH and keeps only putative crossover (CO)-associated LOH
		python ./PhaseClean_RemoveNCOLengthPBs.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' $NCOcutoff
		## create a file that removes all CO-associated LOH and keeps only putative NCO-associated LOH
		python ./PhaseClean_RemoveCOLengthPBs.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' $NCOcutoff
		###### then for the haploids
		## create a file that removes all NCO-associated LOH and keeps only putative crossover (CO)-associated LOH
		python ./PhaseClean_RemoveNCOLengthPBs.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' $NCOcutoff
		## create a file that removes all CO-associated LOH and keeps only putative NCO-associated LOH
		python ./PhaseClean_RemoveCOLengthPBs.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' $NCOcutoff
		## Once you've removed NCOs, there will be some CO-associated haplotype blocks that were previously interrupted by NCOs. These need to be merged.
		python ./PhaseClean_MergeAdjacents.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP.txt'
		## Once you've done that merging, we can tally the number of recombination events detected.
		python ./PhaseClean_CountPhaseSwitches.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents.txt'
		## create a file for HaplotypeBlock presentation
		python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan'$NCOcutoff'BP.txt'
		python ./PhaseClean_RemoveFirstBlockForEachContig.py ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents.txt'
		## create a file that contains all recombination events, including crossovers and gene conversions (NCOs)
		cat ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt' > ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
		echo '' >> ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
		tail -n +2 ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan'$NCOcutoff'BP_FirstBlocksOnContigRemoved.txt' >> ./data/$LINE'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan'$NCOcutoff'BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
	done
done


######################################################################################################
###
###
###  Now we will perform PAIRWISE analyses of recombination across haploid males from multiple colonies of each clonal line
###
###
######################################################################################################


####### Pairwise analyses
# mkdir -p ./data/Pairwise/LineA/DipFemLOH/
# mkdir -p ./data/Pairwise/LineA/HapMaleRec/
# mkdir -p ./data/Pairwise/LineB/DipFemLOH/
# mkdir -p ./data/Pairwise/LineB/HapMaleRec/
# 
# for SAMPLE in SM76 SM83 SM87 KL46 KL66
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM65 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents_ForPairwise.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/Pairwise/LineA/HapMaleRec/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt'
# done
# 
# 
# for SAMPLE in SM74 SM75 SM56 KL70
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM03 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents_ForPairwise.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/Pairwise/LineB/HapMaleRec/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt'
# done
# 
# 
# for SAMPLE in W07 C16B-3-1 C16B-4-2 C18 C1A-1
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name W06 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table_Contigs' ./data/Pairwise/LineA/DipFemLOH/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineA/DipFemLOH/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# done
# 
# for SAMPLE in W05 STC6L STC1_W STC12_W
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name W04 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table_Contigs' ./data/Pairwise/LineB/DipFemLOH/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineB/DipFemLOH/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# done
# 
