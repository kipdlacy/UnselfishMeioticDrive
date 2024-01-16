this directory documents analysis of loss of heterozygosity between pairs of ants from the recently derived transgenic marker lineage
most of these pairs are of known pedigree (i.e., are mother-daughter or sister-sister pairs).
some of the pairs are of unknown pedigree (i.e., comparisons between individiuals from different known pedigrees).
The analysis runs from two different "Controller" bash scripts - one for the short read data, and one for the linked read data.
These "controller" scripts call the GATK, vcftools, and several custom python scripts, which can call be found in this directory.
Note that the "diffs" output files, which contain all of the sites that differ in heterozygosity between pairs of samples
contain false positives, and need to be manually screened (by viewing in IGV) before any particular difference can be trusted.

The "Controller" scripts also contain code to assess the proportion of the genome that contains heterozygous SNPs.
This is achieved by taking the well-filtered heterozygous SNP list for each sample and counting the number of 
10kb nonoverlappint (adjacent) windows that contain at least one heterozygous SNP.
