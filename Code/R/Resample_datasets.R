#===================================================================================================================================
# RESAMPLE SECTION
#===================================================================================================================================

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
