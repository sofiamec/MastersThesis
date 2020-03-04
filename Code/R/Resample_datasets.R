#===================================================================================================================================
# RESAMPLE SECTION
#===================================================================================================================================

# Required from run_entire_analysis script

# Data = Gut2
# saveName = "Gut2"
# saveExpDesign = "m60_d2e6_q10_f10"
# m = 60 nr     samples in one group
# d = 2000000   seq. depth
# q = 10        the fold change for downsampling
# f = 0.1       the fraction of genes to be downsampled

#################################################################################################################

# Resample data
ResampData=resample(Data=Data, m=m, d=d)

# check number of genes wiht low counts in the prduced datasets
countsResampData=compute_low_counts(ResampData)

#################################################################################################################

# Downsampling the resampled dataset
resultList<- introducing_DAGs(Data = ResampData, q = q, f = f)
DownSampledData<-resultList[[1]]
DAGs<-resultList[[2]]

## Strata
if (runStrata==T){
  # Extract and sort gene-strata acording to DAGs
  corrStrata<-DataStrata[rownames(DataStrata) %in% rownames(DAGs),]
  corrStrata<-corrStrata[rownames(DAGs),]
  DAGs<-data.frame(DAGs,corrStrata)
  rm(corrStrata)
}



# Saving downsampled datasets and corresponding overview of DAGs
write.csv(DownSampledData, file=sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run))
write.csv(DAGs, file=sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run))

# remove variables/datasets
rm(ResampData, countsResampData, resultList)
