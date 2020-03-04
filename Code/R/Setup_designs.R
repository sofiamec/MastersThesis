#===================================================================================================================================
# SET-UP SECTION
#===================================================================================================================================

# Setup script for assigning correct names for a given design as well as preprocessing of datasets
# This script has to be run in combination with the "Run_entire_analysis"-script

 
#===================================================================================================================================
# FUNCTIONS SECTION
#===================================================================================================================================

## FUNCTION to Remove low counts 
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

## FUNCTION for creating strata
# DESeq2-analysis for original data 
# This function uses DESeq2 to estimate the mean count and dispersion for genes in a dataset 
# (after normalising them based on sequencing depth). This is then used for creating strata.
#
# Input:  Data = the data to analyse (Gut2 or Marine)
#         numberOfStrata = the number of strata/levels used when dividing the data
# Output: Result = a dataframe containing the mean values and the estimated dispersion for each gene 
#                  as well as the corresponding abundance- and variability-strata for each gene
DESeq2_for_strata=function(Data,numberOfStrata){
  
  DesignMatrix <- data.frame(group=factor(rep(1,(ncol(Data)))))                        # all samples belong to the same group 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data,DesignMatrix, design=~1) # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                           # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq,independentFiltering=FALSE,cooksCutoff=FALSE)         # extract results
  
  Result <- data.frame(rownames(Data), Res$baseMean, dispersions(ResultDESeq))  # dataframe with genes, their mean values (after normalization) and dispersion estimates
  Result <- data.frame(Result[,-1], row.names=Result[,1])                       # put first column (genes) as rowname
  
  # create a vector for assigning smallest to highest strata
  n=numberOfStrata
  repVector<-as.integer(cumsum(c(rep(nrow(Data)/n,n))))
  repVector<-repVector-c(0,repVector)[-(length(repVector)+1)]
  strataVector<-rep(1:n,repVector)
  
  # sort by abundance and assign abundance-strata
  Result<-Result[order(Result[,1]),] 
  Result<-cbind(Result,strataVector)
  
  # sort by variability and assign variability-strata
  Result<-Result[order(Result[,2]),]
  Result<-cbind(Result,strataVector)
  
  # rename columns, factorise strata and sort genes by original order
  colnames(Result)<-c("BaseMean", "Estimated dispersion", "AbundanceStrata", "VariabilityStrata")
  Result$AbundanceStrata<-as.factor(Result$AbundanceStrata)
  Result$VariabilityStrata<-as.factor(Result$VariabilityStrata)
  Result<-Result[rownames(Data),]
  
  rm(n, numberOfStrata, repVector,strataVector)
  return(Result)
}

#===================================================================================================================================
# SET-UP SECTION
#===================================================================================================================================
## Loading original data
#  filtering out samples not meeting expDesign requirements
#  filtering out genes with too low counts

if(saveName == "Gut2"){
  plotName = "Human Gut II"
  Gut2Original <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
  Gut2Intermediate = Gut2Original[,colSums(Gut2Original)>=5000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Gut2 <- remove_low_counts(Gut2Intermediate)
  Data = Gut2
  
  boldvalue2="0"
  relations<-c(3000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000)
    sequencingDepthName<-c("10k","100k","500k","1M", "5M")
  }
  
  rm(Gut2, Gut2Original, Gut2Intermediate)        # remove original and intermediate datasets
  
} else if(saveName == "Marine"){
  plotName = "Marine"
  MarineOriginal <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)
  MarineIntermediate = MarineOriginal[,colSums(MarineOriginal)>=10000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Marine <- remove_low_counts(MarineIntermediate)
  Data = Marine
  
  boldvalue2="5e+07"                                              
  relations<-c(3000000,5000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000,10000000)
    sequencingDepthName<-c("10k","100k","500k","1M", "5M", "10M")
  }
  
  rm(Marine, MarineOriginal, MarineIntermediate)  # remove original and intermediate datasets
}

if (runStrata==T){
  DataStrata<-DESeq2_for_strata(Data,numberOfStrata)
}
#===================================================================================================================================
## Setting upp the right environment and naming designs

AllSaveDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 
AllPlotDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 

for (effect in 1:length(effectsizes)) {           # looping over q
  q=effectsizes[effect]
  
  # Creating the standard designs
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      
      AllSaveDesigns[effect,group,seq] <- sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
      AllPlotDesigns[effect,group,seq] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
      
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group,seq]))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group,seq])), recursive = T)
      }
      
      if (!dir.exists(sprintf("../../Result/%s/IntermediatePlots", saveName))){
        dir.create(file.path("../../Result", sprintf("%s", saveName), "/IntermediatePlots"), recursive = T)
      }
    }
  }
  
  # Creating 3 extra designs                               
  if (extraDesigns==T){
    for (i in 1:extraL) {              # looping over extra designs with fixed m and d
      m=extraGroups[i]
      d=extraSeqDepth[i]
      dD=extraSeqDepthName[i]
      
      AllSaveDesigns[effect,group+i,seq+i] <- sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
      AllPlotDesigns[effect,group+i,seq+i] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
      
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group+i,seq+i]))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group+i,seq+i])), recursive = T)
      }
      
      if (!dir.exists(sprintf("../../Result/%s/IntermediatePlots", saveName))){
        dir.create(file.path("../../Result", sprintf("%s", saveName), "/IntermediatePlots"), recursive = T)
      }
    }}
}

rm(d,dD,effect, group, m, q, seq, remove_low_counts, DESeq2_for_strata)
