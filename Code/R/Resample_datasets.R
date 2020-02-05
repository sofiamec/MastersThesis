# Resample new datasets

# Original datasets
Gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
Gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
Marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)

# In order to get the same results each time
seed=100
set.seed(seed=seed)

#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================
# Resample
# This function resamples a new datasets with n*2 samples and d reads in each sample (depth) 
# input arguments:
  # Data = the data to sample from
  # n = number of samples in the new datasets
  # d = sequencing depth for each sample in the new datasets
# output:  the resampled data in a large dataframe containing n*2 groups
resample = function(Data, n, d){
  sampleVector=sample(ncol(Data), 2*n)                      # vector w. the columnnumber of the sampled samples for both datasets
  DataNew=data.frame(row.names(Data), stringsAsFactors = F) # dataframe to put the resampled data in
  
  for (i in 1:(2*n)){
    readList <- rep(row.names(Data), times=Data[,sampleVector[i]])                # vector containing each read as one entry 
    sampledReads <- sample(readList, size=d)                                      # vector with d resampled reads from "readList"
    sampledVector <- as.data.frame(table(sampledReads), stringsAsFactors = FALSE) # assemble vector "sampledReads" into dataframe
    colnames(sampledVector) <- c("Gene",colnames(Data[sampleVector[i]]))          # name the columns according to the sample-name  
    DataNew <- merge(DataNew, sampledVector, by.x = 1, by.y = 1, all.x = T)       # insert "sampledVector" to "DataNew"
  }
  
  rownames(DataNew) <- DataNew[,1]                  # add genes as row names
  DataNew <- DataNew[,-1]                           # remove column containg genes
  DataNew[is.na(DataNew)] <- 0                      # set all "NA" to 0
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
  cat(sprintf("Number of genes with low counts:        %s\n", sum(r!=0)))
  cat(sprintf("Number of genes wiht acceptable counts: %s\n", sum(r==0)))
  return(r)
}

#===================================================================================================================================

# Remove low counts
# For a given dataset, this function removes genes with low counts (>75 % or an average count <3).
# input: Data = the data to remove genes from
# output: 
remove_low_counts=function(Data){
  
  a=rowSums(Data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(Data)) {
    b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
    r[i]<-a[i]+b[i]
  }
  FilteredData=Data[r==0,]
  return(FilteredData)
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
  
  
  # if we allow unbalanced DAGs (May not work)
  #rG1 <- randomGenes[1:floor(nDAGs/2)] # will be downsampled in dataset 1 
  #rG2 <- randomGenes[ceiling(nDAGs/2):nDAGs] # will be downsampled in dataset 2
  
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

################################### Test whether to filter the original data or not  ##############################
# Remove when desicion has been made!

# WITHOUT filtering the original data
ResampData=resample(Gut2, n=60, d=2000000)  # resampling of data

# check number of genes wiht low counts in the prduced dataset
countsResampData=compute_low_counts(ResampData)  

# filter resampled datasets
ResampDataFilter=remove_low_counts(ResampData)
summary(colSums(ResampDataFilter))

# WITH filtering of the original data
Gut2Filter <- remove_low_counts(Gut2)

ResampData2=resample(Gut2Filter, n=60, d=2000000)  # resampling of data

# check number of genes wiht low counts in the prduced datasets
countsResampData2=compute_low_counts(ResampData2)

# Filter the resampled datasets
ResampData2Filter=remove_low_counts(ResampData2)
summary(colSums(ResampData2Filter))

#################################################################################################################

# Save generated dataset to intermediate folder
write.csv(ResampData2, file=sprintf("../../Intermediate/ResampData_seed%d.csv",seed))

# Load generated dataset from intermediate folder
ResampData2 <- read.csv(file=sprintf("../../Intermediate/ResampData_seed%d.csv",seed), header = T, row.names = 1)

#################################################################################################################

# Downsampling the resampled dataset
resultList<- introducing_DAGs(Data = ResampData2, q = 10, f = 0.10)
downSampledData<-resultList[[1]]
DAGs<-resultList[[2]]

# Saving downsampled datasets and corresponding overview of DAGs
write.csv(downSampledData, file=sprintf("../../Intermediate/downSampledData_seed%d.csv",seed))
write.csv(DAGs, file=sprintf("../../Intermediate/DAGs_seed%d.csv",seed))