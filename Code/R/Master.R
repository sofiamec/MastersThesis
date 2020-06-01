## Selecting parameters and data:
onTerra = T                                                 # use T if running analysis on Terra (large scale settings applied)
saveName = "Resistance"  # "Gut2", "Marine" or "Resistance      # this will in turn load the correct data
f = 0.10                                                    # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
runStrata = F
extraDesigns = T                                              # use T if the analysis of DAGs should be performed with DESeq2. Use F to choose OGLM instead
analyses <- c("DESeq")#, "OGLM")

source("Functions_and_libraries.R")

#########################################


suppressWarnings(source("Setup_designs_final.R"))
  
loadData = T # CHANGE THIS TO "T" IN THE FINAL VERSION!!  
for (analysis in analyses){
  cat(sprintf("==================== %s_%s ====================\n", saveName, analysis), file = stdout()) 
  source("Run_entire_analysis_final.R")
  loadData = T  
}
  
rm(relations, boldvalue2, AllPlotDesigns, plotExpDesign, saveExpDesign, repeats, extraL)
