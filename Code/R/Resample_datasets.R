# Resample new datasets

# Original datasets
Gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
Gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
Marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)


#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================
## Resample-function
  # input arguments:
    # Data = the data to sample from
    # n = number of samples in the new datasets
    # m = sequencing depth for each sample in the new datasets
  # output:  the resampled data in a large dataframe containing n*2 groups
resample = function(Data, n, m){
  sampleVector=sample(ncol(Data), 2*n)                      # vector w. the columnnumber of the sampled samples for both datasets
  DataNew=data.frame(row.names(Data), stringsAsFactors = F) # dataframe to put the resampled data in
  
  for (i in 1:(2*n)){
    readList <- rep(row.names(Data), times=Data[,sampleVector[i]])                # vector containing each read as one entry 
    sampledReads <- sample(readList, size=m)                                      # vector with m resampled reads from "readList"
    sampledVector <- as.data.frame(table(sampledReads), stringsAsFactors = FALSE) # assemble vector "sampledReads" into dataframe
    colnames(sampledVector) <- c("Gene",colnames(Data[sampleVector[i]]))          # name the columns according to the sample-name  
    DataNew <- merge(DataNew, sampledVector, by.x = 1, by.y = 1, all.x = T)       # insert "sampledVector" to "DataNew"
  }
  
  rownames(DataNew) <- DataNew[,1]                  # add genes as row names
  DataNew <- DataNew[,-1]                           # remove column containg genes
  DataNew[is.na(DataNew)] <- 0                      # set all "NA" to 0
  return(list(DataNew[1:n], DataNew[(n+1):(2*n)]))  # didivde into two datasets
}

#===================================================================================================================================

# Summary of low counts
# For a given dataset, this function computes genes which will be removed when filtering.
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
  cat(sprintf("Removed: %s\n", sum(r!=0)))
  cat(sprintf("Kept: %s\n", sum(r==0)))
  return(r)
}

#===================================================================================================================================
#===================================================================================================================================


ResampledData=resample(Gut2, n=50, m=2000000)  # resampling of data
Data1=ResampledData[[1]]                       # extract resampled dataset1
Data2=ResampledData[[2]]                       # extract resampled dataset2

# check number of genes wiht low counts in the prduced datasets
countsData1=compute_low_counts(Data1)
countsData2=compute_low_counts(Data2)

