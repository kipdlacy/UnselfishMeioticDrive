# The purpose of this script is to perform multiple sequence alignments of junction reads
# To do this we are using the R package "msa" (https://bioconductor.org/packages/release/bioc/html/msa.html)

# Installing the package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")

# Loading the package
library(msa)

# Loading the junction reads from the beginning of the insert
B_ie1dsRed_Front <- readDNAStringSet('./data/ie1-dsRed_STC6_ReadsAlignedToInsert_Front.txt')
B_ie1dsRed_Front
# Aligning with CLUSTAL
B_ie1dsRed_Front_Align <- msa(B_ie1dsRed_Front)
B_ie1dsRed_Front_Align
# Printing the multiple sequence alignment
print(B_ie1dsRed_Front_Align, show='complete')

# Loading the junction reads from the end of the insert
B_ie1dsRed_Rear <- readDNAStringSet('../data/ie1-dsRed_STC6_ReadsAlignedToInsert_Rear.txt')
B_ie1dsRed_Rear
# Aligning with CLUSTAL
B_ie1dsRed_Rear_Align <- msa(B_ie1dsRed_Rear)
B_ie1dsRed_Rear_Align
# Printing the multiple sequence alignment
print(B_ie1dsRed_Rear_Align, show='complete')
