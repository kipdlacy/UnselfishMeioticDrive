# load required package
library(ggplot2)

# Loop through each sample
for (SAMPLE in c("MD-B-1-M")) {
  # Read the data from the file into a DataFrame
  # The file name is constructed using the SAMPLE variable
  DataFile <- paste0('./data/', SAMPLE, '_Insert_RawDepth_InsertOnly.txt_DivBy41.0_TGnWGdepthForR.txt')
  Data <- read.table(DataFile, header = TRUE)
  head(Data)
  
  # Set up the PDF device to save the plot (the output file name includes the SAMPLE name)
  PDF_Name <- paste0('./results/DepthPlots/', SAMPLE, '.pdf')
  pdf(PDF_Name, width = 5, height = 5)
  Plot<-ggplot(Data,aes(y=NormalizedDepth, x=EndoOrInsert)) +
    geom_violin()  +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 15) +
    ggtitle(SAMPLE) +
    scale_x_discrete(name="Genomic Region", labels=c("TrustedDiploidSites" = "Diploid Sites","TransgeneInsert" = "Transgene Insert")) +
    scale_y_continuous(name='Normalized Read Depth') 
  print(Plot)
  dev.off()
}
  
##### if your transgene uses an endogenous promoter you can make a plot like this.
# library(ggplot2)
# library(stringr)
# gCamp1<-read.table('../../../../../../Volumes/antqueen/genomics/experiments/analyses/KDL20210903_RelativeDepthTransgeneInsertions/results/gCamp1_Insert_RawDepth_InsertOnly.txt_DivBy44.0_ORCO_Insert_RefGenome_Sites4RwithRemovedIQROutliers_Types.txt',header = T)
# head(gCamp1)
# gCamp1
# pdf('./results/DepthPlots/gCamp1.pdf',width = 8,height=5)
# Plot<-ggplot(gCamp1,aes(y=NormalizedDepth, x=EndoOrInsert, color=AlignedTo)) +
#   geom_violin()  +
#   geom_boxplot(width=0.1) +
#   theme_bw(base_size = 15) +
#   ggtitle("Orco>gCamp") +
#   scale_x_discrete(name="Genomic Region", labels=c("EndogenousGenome" = "Heterozygous SNPs","TransgeneInsert" = "Non-ORCO Insert", "ORCOpromoter" = "ORCO Promoter")) +
#   scale_y_continuous(name='Normalized Read Depth') +
#   scale_color_brewer(palette="Set1")
# Plot
# dev.off()