# Resample new datasets

# Original datasets
gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)


#===================================================================================================================================
#=========================================== Function ==============================================================================
#===================================================================================================================================
## Resample-function
  # input arguments:
    # data = the data to sample from
    # n = number of samples in the new datasets
    # m = sequencing depth for each sample in the new datasets
  # output:  the resampled data in a large dataframe containing n*2 groups
resample = function(data, n, m){
  samplevector=sample(ncol(data), 2*n)                      # vector w. the columnnumber of the sampled samples for both datasets
  datanew=data.frame(row.names(data), stringsAsFactors = F) # dataframe to put the resampled data in
  
  for (i in 1:(2*n)){
    readlist <- rep(row.names(data), times=data[,samplevector[i]])                # vector containing each read as one entry 
    sampledreads <- sample(readlist, size=m)                                      # vector with m resampled reads from "readlist"
    sampledvector <- as.data.frame(table(sampledreads), stringsAsFactors = FALSE) # assemble vector "sampledreads" into dataframe
    colnames(sampledvector) <- c("Gene",colnames(data[samplevector[i]]))          # name the columns according to the sample-name  
    datanew <- merge(datanew, sampledvector, by.x = 1, by.y = 1, all.x = T)       # insert "sampledvector" to "datanew"
  }
  
  rownames(datanew) <- datanew[,1]                  # add genes as row names
  datanew <- datanew[,-1]                           # remove column containg genes
  datanew[is.na(datanew)] <- 0                      # set all "NA" to 0
  return(list(datanew[1:n], datanew[(n+1):(2*n)]))  # didivde into two datasets
}

#===================================================================================================================================

resampleddata=resample(marine, n=50, m=2000000)  # resampling of gut1

data1=resampleddata[[1]]      # extract resampled dataset1
data2=resampleddata[[2]]      # extract resampled dataset2

colSums(data1)
