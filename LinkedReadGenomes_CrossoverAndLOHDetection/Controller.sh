#! /usr/bin/bash

### First remove sites with indications of poor genome assembly from putative crossover event lists
### And also create lists of sites with such indications for later screening
bash ./Clean_RemoveBadSitesFromFiles.sh

### Then compare phasing between pairs of individuals and screen such sites and compile all putative crossover events for manual screening
bash ./Compare_PhasingBetweenIndividuals.sh

### Then manually screen all putative crossovers.

### Then obtain the phase of all samples at high quality crossovers
bash ./ObtainPhaseOfAllSamplesAtAllDetectedCrossovers.sh

### Then get the hypothetical proportion of the genome at which we ought to be able to detect crossovers using our method
bash ./DetectableGenomePropGetterController.sh

### Plot data
mkdir -p ./results/Plots/
cd ./results/Plots
R ../../CalcStatsAndPlotLinkedReadCrossoverRecomb.R