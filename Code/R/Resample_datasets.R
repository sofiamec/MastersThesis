# Required from run_entire_analysis script

# Data = Gut2
# saveName = "Gut2"
# saveExpDesign = "m60_d2e6_q10_f10"
# m = 60 nr     samples in one group
# d = 2000000   seq. depth
# q = 10        the fold change for downsampling
# f = 0.1       the fraction of genes to be downsampled


# In order to get the same results each time
rep=run 

#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================
# Resample
# This function resamples a new datasets with m*2 samples and d reads in each sample (depth) 
# input arguments:
  # Data = the data to sample from
  # m = number of samples in the new datasets
  # d = sequencing depth for each sample in the new datasets
# output:  the resampled data in a large dataframe containing m*2 groups
resample = function(Data, m, d){
  sampleVector=sample(ncol(Data), 2*m)                      # vector w. the columnnumber of the sampled samples for both datasets
  DataNew=data.frame(row.names(Data), stringsAsFactors = F) # dataframe to put the resampled data in
  
  for (i in 1:(2*m)){
    readList <- rep(row.names(Data), times=Data[,sampleVector[i]])                # vector containing each read as one entry 
    sampledReads <- sample(readList, size=d)                                      # vector with d resampled reads from "readList"
    sampledVector <- as.data.frame(table(sampledReads), stringsAsFactors = FALSE) # assemble vector "sampledReads" into dataframe
    colnames(sampledVector) <- c("Gene",colnames(Data[sampleVector[i]]))          # name the columns according to the sample-name  
    DataNew <- merge(DataNew, sampledVector, by.x = 1, by.y = 1, all.x = T)       # insert "sampledVector" to "DataNew"
  }
  
  rownames(DataNew) <- DataNew[,1]                  # add genes as row names
  DataNew <- DataNew[,-1]                           # remove column containg genes
  DataNew[is.na(DataNew)] <- as.integer(0)          # set all "NA" to 0 (as integers, since DESeq2 require integer counts)  return(DataNew)
  
  return(DataNew) 
  }

#===================================================================================================================================

# Summary of low counts
# For a given dataset, this function computes genes which should be removed when filtering.
# It also prints how many genes will be removed respectively kept.
# input: Data = the data to be filtered
# output: r = a logical vector with 1 if a gene should be removed, and 0 if it should be kept. 
compute_low_counts=function(Data){
  
  a=rowSums(Data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(Data)) {
    b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
    r[i]<-a[i]+b[i]
  }
  cat(sprintf("Number of genes with low counts:                 %s\n", sum(r!=0)))
  cat(sprintf("Number of genes with acceptable counts:          %s\n", sum(r==0)))
  return(r)
}

#===================================================================================================================================

# Introducing DAGs
# For a resampled dataset (including both groups), this function introduces DAGs by
# downsampling a given fraction of genes. The DAGs will be balanced in the two groups.
# Input: Data = the resampled dataset including both groups (dataset 1 and 2)
#        q = the fold change. Ex. 10. Will result in relative abundance between datasets.
#        f = the desired (total) fraction of genes to be downsampled (the output will not follow this exactly)
# Outut: downSampledData = the new dataset with introduced DAGs. Analysis should be performed on this dataset
#        DAGs = a matrix with an overview of which genes have been downsampled in which dataset with the given q
introducing_DAGs = function(Data, q, f){
  
  downSampledData = Data
  nDAGs = 2 * round(f*nrow(Data)/2) # the total number of genes to be downsampled if we don't allow unbalanced DAGs
  #nDAGs = trunc(f*nrow(Data)+0.5) # the total number of genes to be downsampled if we allow unbalanced DAGs
  
  # Creating empty matrix for overview of DAGs
  DAGs = matrix(ncol = 2, nrow = nrow(Data)) 
  colnames(DAGs) = c(sprintf("Sample 1 to %d", ncol(Data)/2),sprintf("Sample %d to %d", ncol(Data)/2+1, ncol(Data)))
  rownames(DAGs) <- rownames(Data)
  # Selecting random genes
  randomGenes <- sample(nrow(Data),nDAGs) # Selects n random genes in the dataset which will be downsampled. 
  if (nDAGs==1|| nDAGs==0) {
    print("Too few DAGs introduced")
    nDAGs=0
  }
  # if we don't allow unbalanced DAGs
  rG1 <- randomGenes[1:(nDAGs/2)] # will be downsampled in dataset 1 
  rG2 <- randomGenes[(nDAGs/2+1):nDAGs] # will be downsampled in dataset 2
  
  for (gene in rG1) {
    for (sample in 1:(ncol(Data)/2)) {
      downSampledData[gene,sample] <- rbinom(n = 1 ,size = Data[gene,sample] ,prob = 1/q)
    }
    DAGs[gene, 1] <- -q
  }
  
  for (gene in rG2) {
    for (sample in (ncol(Data)/2+1):(ncol(Data))) {
      downSampledData[gene,sample] <- rbinom(n = 1 ,size = Data[gene,sample] ,prob = 1/q) 
      DAGs[gene, 2] <- -q
    }
  }
  
  DAGs<-DAGs[rowSums(DAGs, na.rm=T)!=0,]
  return(list(downSampledData,DAGs))
}

#===================================================================================================================================
#===================================================================================================================================

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

# Saving downsampled datasets and corresponding overview of DAGs
write.csv(DownSampledData, file=sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, rep))
write.csv(DAGs, file=sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, rep))

# remove variables/datasets
rm(ResampData, countsResampData, resultList)
