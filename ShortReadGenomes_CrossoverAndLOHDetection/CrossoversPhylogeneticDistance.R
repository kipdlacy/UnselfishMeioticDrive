setwd("/Users/kiplacy/../../../../Volumes/antqueen/genomics/experiments/analyses/KDL20220511_HapMaleRecDipFemLOH/")

### To determine whether the abundance of crossovers observed among haploid male genomes have accumulated over many generations
### we investigated the association of the number of crossovers with the phylogenetic distance between samples
### To ensure independence of these data, we investigated samples from multiple field-collected (and somewhat genetically distinct)
### colonies from each clonal line, and calculated phylogenetic distance among the diploid females.
### We then analyzed how that scaled with the normalized number of crossovers.

## Loading the following R libraries
library(gdsfmt)
library(SNPRelate) 
library(ggplot2)
library(ggtree)
library(ape)

### Save the path to the GATK vcf
LineADips<-"./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf"
### turn the VCF file into a less data intensive form (GDS) for easier computing
snpgdsVCF2GDS(LineADips,"./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.gds",method ="biallelic.only")
### Prepare the data so it is formatted correctly to create a dissimilarity matrix.
genofileLineADips <- snpgdsOpen("./data/LineAmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.gds")
set.seed(100) ## making the code reproducible
## Test Begin
snpgdsDiss(genofileLineADips,num.thread=2, autosome.only=FALSE)
snpgdsIBS(genofileLineADips,num.thread=2, autosome.only=FALSE)
## Test End
ibs_LineADips <- snpgdsHCluster(snpgdsIBS(genofileLineADips,num.thread=2, autosome.only=FALSE))
rvLineADips <- snpgdsCutTree(ibs_LineADips)
### Saving the dendrograms to new Variables
treeLineADips = rvLineADips$dendrogram
plot(rvLineADips$dendrogram,horiz=T, main ="DipFems, Line A" )
### converting dendrograms to class hclust and to new variables
hcLineADips = as.hclust(rvLineADips$dendrogram)
### Making the hclust object into a phylo object in ape
LineADipsPhylo <- as.phylo(hcLineADips) 
#### Writing to a newick tree file
write.tree(phy=LineADipsPhylo , file="./results/LineADipFemsPhylo.newick") 
### DIstances
dist.nodes(LineADipsPhylo)
cophenetic.phylo(LineADipsPhylo)


### Save the path to the GATK vcf as a variable
LineBDips<-"./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf"
### turn the VCF file into a less data intensive form (GDS) for easier computing
snpgdsVCF2GDS(LineBDips,"./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.gds",method ="biallelic.only")
### Prepare the data so it is formatted correctly to create a dissimilarity matrix.
genofileLineBDips <- snpgdsOpen("./data/LineBmultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.gds")
set.seed(100) ## making the code reproducible
ibs_LineBDips <- snpgdsHCluster(snpgdsIBS(genofileLineBDips,num.thread=2, autosome.only=FALSE))
rvLineBDips <- snpgdsCutTree(ibs_LineBDips)
### Saving the dendrograms to new Variables
treeLineBDips = rvLineBDips$dendrogram
plot(rvLineBDips$dendrogram,horiz=T, main ="DipFems, Line B" )
### converting dendrograms to class hclust and to new variables
hcLineBDips = as.hclust(rvLineBDips$dendrogram)
### Making the hclust object into a phylo object in ape
LineBDipsPhylo <- as.phylo(hcLineBDips) 
#### Writing to a newick tree file
write.tree(phy=LineBDipsPhylo , file="./results/LineBDipFemsPhylo.newick") 
### DIstances
dist.nodes(LineBDipsPhylo)
cophenetic.phylo(LineBDipsPhylo)


######################################################################################################################
###################       CoPlotting the data from both clonal lines using results from above          ###############
######################################################################################################################

### First we obtain genetic distances

# These are the phylogenetic distances among STC6 females
DistAmongSTC6DipFems=c(0.002764319,0.002764319,0.0003313041)
# Find the mean among STC6 females and use that as the pairwise distance among STC6 samples
MeanDistAmongSTC6DipFems=mean(DistAmongSTC6DipFems)
# Record the phylogenetic distances between STC1 and STC6, and between STC12 and STC6
DistBtwnSTC1nSTC6=0.05565909
DistBtwnSTC12nSTC6=0.05565909
# These are the phylogenetic distances among C16 females 
DistAmongC16DipFems=c(0.02840124,0.02840124,0.02840124,0.0065823550,0.0005611188,0.006582355)
# Find the mean
MeanDistAmongC16DipFems=mean(DistAmongC16DipFems)
# Record the phylogenetic distances between C1 and C16, and between C18 and C16
DistBtwnC18nC16=0.08973584
DistBtwnC1nC16=0.10227901
# Record the number of crossovers recorded between each sample and the focal sample for each colony
# For C16, the focal individual is SM65. For STC6, the focal individual is SM03
NumCOs=c(157,302,1790,1827,171,172,181,595,696)
# Because crossover detectability using comparisons among whole genome sequences scales with the number of genetic markers
# we need to normalize by the number of ancestrally heterozygous SNPs in each clonal line and colony
NumSNPsLineAHap=11160
NumSNPsLineBHap=47017
# After dividing the number of crossovers by the Number of SNPs for each clonal line
NumCOProps=c(0.003339218,0.006423209,0.038071336,0.038858285,0.015322581,0.015412186,0.016218638,0.053315412,0.062365591)
# Record the sample names
SampNames=c('SM74','SM75','KL70','SM56','SM76','SM83','SM87','KL46','KL66')
# Record the colonies (with letters appended that are useful for plotting)
COLs=c('D_STC6','D_STC6','F_STC1','E_STC12','A_C16','A_C16','A_C16','C_C1','B_C18')
# Col Type, where 0 means the sample is from the same colony as the focal sample, and 1 and 2 are the other colonies
ColType=c(0,0,1,2,0,0,0,1,2)
# The clonal line of each colony
CloneLine=c("B","B","B","B","A","A","A","A","A")
# Assemble the distances into a vector
Dists=c(MeanDistAmongSTC6DipFems,MeanDistAmongSTC6DipFems,DistBtwnSTC1nSTC6,DistBtwnSTC12nSTC6,MeanDistAmongC16DipFems,MeanDistAmongC16DipFems,MeanDistAmongC16DipFems,DistBtwnC1nC16,DistBtwnC18nC16)
# bind into a dataframe
data=cbind(SampNames,COLs,ColType,Dists,NumCOs,CloneLine,NumCOProps)
data=as.data.frame(data)

### Plotting
pdf("/Users/kiplacy/Desktop/BothLinesHapRec_vs_PhyDist_ColCol_NoSM06_Prop.pdf",width = 5.5,height = 4)
ggplot(data=data,aes(x=as.numeric(Dists),y=as.numeric(NumCOProps),color=CloneLine)) +
  geom_jitter(size=5, shape=ColType,stroke=1.5) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  theme_bw(base_size = 15) +
  xlab('Pairwise Distance with Focal Colony') +
  ylab('Number of Crossovers / Number of SNPs') +
  scale_color_manual(values=c("purple","orange"))
dev.off()

#####################
### Linear Models ###
#####################

## Split data by clonal line
CloneData<-split(data,data$CloneLine)

# Look at the Line A side
CloneData$A
# Inspect the correlation
Acorr=cor(as.numeric(CloneData$A$Dists),as.numeric(CloneData$A$NumCOProps))
Acorr
# Assign phylogenetic distance to the "X" value
A_x<-as.numeric(CloneData$A$Dists)
# Assign normalized crossover number to the "Y" value
A_y<-as.numeric(CloneData$A$NumCOProps)
# Bind to a dataframe
Aframe=cbind(A_x,A_y)
Aframe<-as.data.frame(Aframe)
# Construct a linear model
AlinearModel<-lm(A_y ~ A_x, data=Aframe)
# Inspect the linear model
summary(AlinearModel)

# Look at the Line B side
CloneData$B
# Inspect the correlation
Bcorr=cor(as.numeric(CloneData$B$Dists),as.numeric(CloneData$B$NumCOProps))
Bcorr
# Assign phylogenetic distance to the "X" value
B_x<-as.numeric(CloneData$B$Dists)
# Assign normalized crossovers to the "Y" value
B_y<-as.numeric(CloneData$B$NumCOProps)
# Bind to a dataframe
Bframe=cbind(B_x,B_y)
Bframe<-as.data.frame(Bframe)
# Construct a linear model
BlinearModel<-lm(B_y ~ B_x, data=Bframe)
# Inspect the linear model
summary(BlinearModel)
