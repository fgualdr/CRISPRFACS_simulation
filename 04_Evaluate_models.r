#!/usr/bin/env Rscript 
# qsub -I -l select=1:ncpus=5:mem=20gb
# singularity shell -B /hpcnfs docker://fgualdr/crispr_tools 
# suppress verbosity when loading libraries:
options(warn=-1)
library(PRROC)

# Get all the results from the screen
meta_simulations <- "/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/meta_data.txt"
meta_data <- read.delim(meta_simulations, header = TRUE, sep = "\t")

# the results are within the folder: /hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/Simulated_data_counts/
# within there are subfolders all with name /test_XXX/ where XXX is the number of the simulation 
# then within we have "CrisprReshigh_vs_low/Testing_stat/Multivariate_high_vs_low_prob_per_TARGET.txt" which has the result of the tests

# Aim is to generate the ROC and PR AUC for each simulation and add the computed values to the column in meta_data# then we can plot the global result displaying globally the performances and the impact of the individual parameters



# Then we can test the performance of sets of other tools: 
# CRISPRBetaBinomial
# CRISPhieRmix
# CRISPRcleanR
# gscreend
# MAGeCKFlute
