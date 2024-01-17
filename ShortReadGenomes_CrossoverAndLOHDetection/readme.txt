This directory contains scripts to identify crossover events among diploid and haploid genomes from unknown pedigrees.
For diploids, crossovers can only be inferred indirectly using losses of heterozygosity as a proxy.
To detect those we identify runs of homozygosity.
For haploids, we can identify crossovers directly as changes in relative phasing of sites that are heterozygous in diploids from the same clonal background.
The controller script contains the commands to run this analysis, and several steps utilize custom python scripts.
Plotting is done in R scripts.
