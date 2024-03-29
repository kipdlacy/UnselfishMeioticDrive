setwd("/Users/kiplacy/../../../../Volumes/antqueen/booster/KDL20240109_UploadCodeParthPaper/UnselfishMeioticDrive/ShortReadGenomes_CrossoverAndLOHDetection/")
getwd()

### This script is designed to plot recombination events observed among haploid male genomes, 
### and inferred among diploid female genomes using segmental losses of heterozygosity as a proxy.

###########################################
#                                         #
#                                         #
#     Plotting crossovers genome-wide     #
#                                         #
#                                         #
###########################################

### Load plotting packages
library(ggplot2)
library(karyoploteR)

### Import genome-specific files for KaryoploteR
Obir5.4withChr<-toGRanges("../UsefulFiles/GenomeForKaryoplotR_Obir5_4_withChr.txt")
Obir5.4withChr
ObirTigswithChr<-toGRanges("../UsefulFiles/ContigsForKaryoplotR_Obir5_4_withChr.txt")
ObirTigswithChr

### Import data files
# Crossovers inferred among C16 diploids using LOH as a proxy
C16_LOH_CO<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_NoPBsLessThan1500BP_Curated_COonly.txt')
# Crossovers observed among C16 haploids
C16_HapRec<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_PSforGRanges_ManuallyScreened.txt')
# Putatively ancestrally heterozygous sites in colony C16
C16_PutAncHet<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener')
# Crossovers inferred among STC6 diploids using LOH as a proxy
STC6_LOH_CO<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_NoPBsLessThan1500BP_Curated_COonly.txt')
# Crossovers observed among STC6 haploids
STC6_HapRec<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_PSforGRanges_ManuallyScreened.txt')
# Putatively ancestrally heterozygous sites in colony STC6
STC6_PutAncHet<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener')

### Plot this data genome-wide using Karyoploter
pdf("./results/Karyoplot_BothLines_UnknownPedigree_Border_LOHCOs.pdf",width = 15,height = 5)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
at <- autotrack(current.track = 6,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=C16_LOH_CO,col='black',border='black')
kpAddLabels(kp, labels="LOH among Diploid Females", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 5,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=C16_HapRec,col='black',border='black')
kpAddLabels(kp, labels="CO among Haploid Males", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=C16_PutAncHet,window.size = 200000,col='grey',border=NA)
kpAxis(kp, ymax=65, r0=at$r0, r1=at$r1, cex=0.8)
at <- autotrack(current.track = 3,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=STC6_LOH_CO,col='black',border='black')
kpAddLabels(kp, labels="LOH among Diploid Females", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=STC6_HapRec,col='black',border='black')
kpAddLabels(kp, labels="CO among Haploid Males", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=STC6_PutAncHet,window.size = 200000,col='grey',border=NA)
kpAxis(kp, ymax=151, r0=at$r0, r1=at$r1, cex=0.8)
dev.off()

### Get axes for density plots
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=C16_PutAncHet,window.size = 200000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density ## result: 65
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=STC6_PutAncHet,window.size = 200000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density ## result: 151

#########################################################################
#                                                                       #
#                                                                       #
#     Haplotype block length distributions for Colonies C16 and STC6    #
#                                                                       #
#                                                                       #
#########################################################################

library(ggplot2)
library(scales)
options(scipen=0)

####################### Line A Colony C16

### Load the haplotype blocks for recombination events observed among C16 haploids
C16HapRecHaplotypes<-read.table("./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt",header = T)
C16blockLengths<-as.data.frame(C16HapRecHaplotypes$InfBlockLength)
C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`<-as.numeric(C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)
### Load the LOH blocks for recombination events inferred via LOH among C16 diploids
C16DipFemHaplotypes<-read.table("./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_ManuallyScreened.txt",header = T)
C16DipFemblockLengths<-as.data.frame(C16DipFemHaplotypes$InfLOHLength)
C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`<-as.numeric(C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`)

### Plot haplotype block length histogram
pdf("./results/BlockLengths_C16haprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",limits = c(-100000,4000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### C16 recombination among haploids--Plot haplotype block length vs SNP number scatterplot to show that long haplotype blocks are not supported by just a few SNPs
pdf("./results/NumSNPsVBlockLengths_C16haprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(C16HapRecHaplotypes, aes(x=C16HapRecHaplotypes$InfBlockLength,y=C16HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,1000000,2000000,3000000,4000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Plot LOH block length histogram
pdf("./results/BlockLengths_C16dipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(C16DipFemblockLengths, aes(x=C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)", breaks = c(0,500000,1000000,1500000,2000000,2500000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### C16 LOH among diploids -- scatterplot of block length vs SNP number to show that long LOH tracts are supported by a corresponding number of SNPs
pdf("./results/NumSNPsVBlockLengths_C16dipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(C16DipFemHaplotypes, aes(x=C16DipFemHaplotypes$InfLOHLength,y=C16DipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000,2000000,2500000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()

####################### Line B Colony STC6

### Load the haplotype blocks for recombination events observed among STC6 haploids
STC6HapRecHaplotypes<-read.table("./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions_ManuallyScreened.txt",header = T)
STC6blockLengths<-as.data.frame(STC6HapRecHaplotypes$InfBlockLength)
STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`<-as.numeric(STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)
### Load the LOH blocks for recombination events inferred via LOH among STC6 diploids
STC6DipFemHaplotypes<-read.table("./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_ManuallyScreened.txt",header = T)
STC6DipFemblockLengths<-as.data.frame(STC6DipFemHaplotypes$InfLOHLength)
STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`<-as.numeric(STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`)

### STC6 haploids -- Plot haplotype block length histogram
pdf("./results/BlockLengths_STC6haprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(STC6blockLengths, aes(x=STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",limits = c(-100000,6000000),breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### STC6 recombination among haploids--Plot haplotype block length vs SNP number scatterplot to show that long haplotype blocks are not supported by just a few SNPs
pdf("./results/NumSNPsVBlockLengths_STC6haprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(STC6HapRecHaplotypes, aes(x=STC6HapRecHaplotypes$InfBlockLength,y=STC6HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### STC6 diploids -- Plot LOH block length histogram
pdf("./results/BlockLengths_STC6dipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(STC6DipFemblockLengths, aes(x=STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000), limits = c(-50000,1500000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### STC6 LOH among diploids -- scatterplot of block length vs SNP number to show that long LOH tracts are supported by a corresponding number of SNPs
pdf("./results/NumSNPsVBlockLengths_STC6dipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(STC6DipFemHaplotypes, aes(x=STC6DipFemHaplotypes$InfLOHLength,y=STC6DipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000,2000000), limits=c(-50000,1500000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()

#########################################################################
#                                                                       #
#                                                                       #
#     Haplotype block length distributions for clonal lines A and B     #
#                                                                       #
#                                                                       #
#########################################################################

####################### Multicolony Line A

### Recombination among haploids -- phase blocks from Line A, multi-colony dataset
LineAHapRecHaplotypes<-read.table("./data/LineA_MultiCol_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt",header = T)
LineAblockLengths<-as.data.frame(LineAHapRecHaplotypes$InfBlockLength)
LineAblockLengths$`LineAHapRecHaplotypes$InfBlockLength`<-as.numeric(LineAblockLengths$`LineAHapRecHaplotypes$InfBlockLength`)
### Loss of heterozygosity among diploids -- LOH blocks from Line A, multi-colony dataset
LineADipFemHaplotypes<-read.table("./data/LineA_MultiCol_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt",header = T)
LineADipFemblockLengths<-as.data.frame(LineADipFemHaplotypes$InfLOHLength)
LineADipFemblockLengths$`LineADipFemHaplotypes$InfLOHLength`<-as.numeric(LineADipFemblockLengths$`LineADipFemHaplotypes$InfLOHLength`)

### Histogram of recombination among line A haploid block lengths
pdf("./results/MultiColonyNonPairwise/BlockLengths_LineAhaprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(LineAblockLengths, aes(x=LineAblockLengths$`LineAHapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block",breaks = c(0,250000,500000,750000,1000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### Scatterplot of recombination haplotype block lengths plotted against the number of constituent SNPs
pdf("./results/MultiColonyNonPairwise/NumSNPsVBlockLengths_LineAhaprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(LineAHapRecHaplotypes, aes(x=LineAHapRecHaplotypes$InfBlockLength,y=LineAHapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block",breaks = c(0,250000,500000,750000,1000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Histogram of LOH among line A diploid block lengths
pdf("./results/MultiColonyNonPairwise/NumSNPsVBlockLengths_LineAdipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(LineADipFemHaplotypes, aes(x=LineADipFemHaplotypes$InfLOHLength,y=LineADipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Loss-of-Heterozygosity Tract",breaks = c(0,1000000,2000000,3000000,4000000)) +
  scale_y_continuous(name="Number of SNPs in Tract") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Scatterplot of LOH block lengths plotted against the number of constituent SNPs
pdf("./results/MultiColonyNonPairwise/BlockLengths_LineAdipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(LineADipFemblockLengths, aes(x=LineADipFemblockLengths$`LineADipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name = "Length of Loss-of-Heterozygosity Tract",breaks = c(0,1000000,2000000,3000000,4000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

####################### Multicolony Line B

### Recombination among haploids -- phase blocks from Line A, multi-colony dataset
LineBHapRecHaplotypes<-read.table("./data/LineB_MultiCol_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt",header = T)
LineBblockLengths<-as.data.frame(LineBHapRecHaplotypes$InfBlockLength)
LineBblockLengths$`LineBHapRecHaplotypes$InfBlockLength`<-as.numeric(LineBblockLengths$`LineBHapRecHaplotypes$InfBlockLength`)
### Loss of heterozygosity among diploids -- LOH blocks from Line A, multi-colony dataset
LineBDipFemHaplotypes<-read.table("./data/LineB_MultiCol_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt",header = T)
LineBDipFemblockLengths<-as.data.frame(LineBDipFemHaplotypes$InfLOHLength)
LineBDipFemblockLengths$`LineBDipFemHaplotypes$InfLOHLength`<-as.numeric(LineBDipFemblockLengths$`LineBDipFemHaplotypes$InfLOHLength`)

### Histogram of recombination among line B haploid block lengths
pdf("./results/MultiColonyNonPairwise/BlockLengths_LineBhaprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(LineBblockLengths, aes(x=LineBblockLengths$`LineBHapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block",breaks = c(0,1000000,2000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### Scatterplot of recombination haplotype block lengths plotted against the number of constituent SNPs
pdf("./results/MultiColonyNonPairwise/NumSNPsVBlockLengths_LineBhaprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(LineBHapRecHaplotypes, aes(x=LineBHapRecHaplotypes$InfBlockLength,y=LineBHapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block",breaks = c(0,1000000,2000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Histogram of LOH among line B diploid block lengths
pdf("./results/MultiColonyNonPairwise/NumSNPsVBlockLengths_LineBdipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(LineBDipFemHaplotypes, aes(x=LineBDipFemHaplotypes$InfLOHLength,y=LineBDipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Loss-of-Heterozygosity Tract",breaks = c(0,1000000,2000000,3000000)) +
  scale_y_continuous(name="Number of SNPs in Tract") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Scatterplot of LOH block lengths plotted against the number of constituent SNPs
pdf("./results/MultiColonyNonPairwise/BlockLengths_LineBdipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(LineBDipFemblockLengths, aes(x=LineBDipFemblockLengths$`LineBDipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name = "Length of Loss-of-Heterozygosity Tract",breaks = c(0,1000000,2000000,3000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

