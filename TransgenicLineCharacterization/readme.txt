This directory documents code used to characterize the transgenic marker lineage that expresses dsRed.
"FitnessTestAndPlot.R" and "SourceData_TransgenicLineFitness.txt" test whether the transgenic line has reduced fitness compared with wild types.

"Controller_TransgeneInsertionNumber.sh" documents analyses to prepare data comparing the relative read depths of regions in the transgene insert and diploid sites in the nuclear genome to count the number of transgene insertions.
Those data are then plotted in "TransgeneInsertionNumberViaDepthPlot.R"
"Controller_TransgeneInsertionNumber.sh" calls several custom python scripts that are also present in this directory.

"CommandsToMapTransgeneInsertionSite.sh" describes how we mapped the insertion site, which included performing a multiple sequence alignment in the script "MSA.R"
