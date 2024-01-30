This directory includes analyses of linked-read whole genome sequencing data of known-pedigree pairs of diploid females. 
Initial processing of data was performed using Universal Sequencing Technologies' tellread and tellsort pipelines, using the commands in "RunTELLREADandTELLSORT.sh"
Then, using the scripts called in Controller.sh, I obtained clean VCF files, identified genomic regions with signatures of poor genome assembly (and therefore low phasing quality),
identified differences in phasing between files (putative crossover events), and screened out such events that are likely to be affected by poor phasing quality. 
Once you have manually screened all of these, you can create a table using ObtainPhaseOfAllSamplesAtAllDetectedCrossovers.sh.
You can obtain the "proportion of the genomes of each pair for which crossovers (ought to be) detectable" using "DetectableGenomePropGetterController.sh"
