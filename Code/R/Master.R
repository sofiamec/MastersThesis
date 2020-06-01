# This script combines the relevant scripts for performing the entire analysis of one specified metagenomic dataset.

# In addition to this scrip, the scripts "Exploration_of_data.R" and "Plot_tables_and_cost.R" may be run to produce additional results.
# Exploration_of_data: Explores the sequencing depth, group sizes, abundance and variability of the real metagenomic datasets.
# Plot_tables_and_cost: May be run once this script has been run for each of the three datasets and will the produce tables for the results
#                       as well as example plots and the results from the economic analysis.

####################  ONLY CHANGE SETTINGS IN THIS SECTION  #####################
source("Package_installations.R")                               # Mute after this script has been run once. It will install all required packages

## Selecting parameters and data:
onTerra = T                                                     # Use T if running analysis on Terra (large scale settings will be applied). 
                                                                # If F is used, other settings have to be specified in "Setup_designs.R" 

saveName = "Resistance"  # "Gut2", "Marine" or "Resistance      # Select on eof the 3 datasets. This will in turn load the correct data
runStrata = T                                                   # Use T if an analysis of strata should be performed
extraDesigns = T                                                # Use T if the analysis of DAGs should be performed with DESeq2. Use F to choose OGLM instead
analyses <- c("DESeq", "OGLM")                                  # Specify if both DESeq and OGLM analyses should be performed

f = 0.10                                                        # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced


###############################################################################

# Loads the required functions and libraries 
source("Functions_and_libraries.R")

# Specifies the experimental designs that should be used
suppressWarnings(source("Setup_designs.R"))
  
# Runs the specified analyses. In the first analysis, datasets are created. In the second analysis, the datasets are loaded 
loadData = F 
for (analysis in analyses){
  cat(sprintf("==================== %s_%s ====================\n", saveName, analysis), file = stdout()) 
  source("Analysis.R")
  loadData = T  
}
  
# Removes unnecessary parameters
rm(relations, boldvalue2, AllPlotDesigns, plotExpDesign, saveExpDesign, repeats, extraL)
