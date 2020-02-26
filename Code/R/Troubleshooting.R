# Troubleshooting high sequencing depths 
library(DESeq2)

#=================================== Choose dataset to troubleshoot ===========================================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 50              # Number of samples in each group (total nr samples = 2*m)
d = 10000000        # Desired sequencing depth per sample. It will not be exct
q = 1.5   #2.5          # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=7    #10           # choose wich run to extract


{ # Quickly gives the case the correct name
  if (saveName == "Gut2"){
    plotName = "Human Gut II"
  } else if (saveName == "Marine"){
    plotName = "Marine"
  } else {
    sprintf("Missing name for dataset")
  }
  
  if (d==1e4||d==1e5||d==5e5){
    dD=d/1000
    prefix="k"
  } else if(d==1e6||d==5e6||d==10e6){
    dD=d/1000000
    prefix="M"
  } else {
    dD=d
    sprintf("wrong d")
    prefix=""
  }
  
  # Names for a certain dataset and name    # Results in:
  saveExpDesign = sprintf("m%d_d%d%s_10q%d_f%d", m, dD, prefix, q*10, f*100)
  plotExpDesign = sprintf("m=%d, d=%d%s, q=%g, f=%d%%",m,dD,prefix,q,f*100)
  
  rm(dD,prefix)
}

# Read dataset 
DownSampledData=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

# ================================================================================================================================

# run dataset in Analysis of DAGs
source("Analysis_of_DAGs_separate.R")


#=================================================================================================================================
# DESeq2-analysis
# This function uses DESeq2 to identfy DAGs in a dataset containing two groups
# Input: Data = the data to analyse
# Output: a dataframe containing the p-value for each gene, ordered with increading p-values
DESeq2_analysis=function(Data){
  
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                     # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(Data), Res$padj, Res$log2FoldChange )          # dataframe with genes, their adjusted p-values and log2fold change
  ResultSorted=as.data.frame(Result[order(Result[,2]),])                    # order with increasing p-value
  ResultSorted <- data.frame(ResultSorted[,-1], row.names=ResultSorted[,1]) # put first column (genes) as rowname
  colnames(ResultSorted) <- c("adjusted p-value", "log2 fold change")       # name the columns
  
  return(ResultSorted)
}

#=================================================================================================================================
# check summary for TP and FP
Result=DESeq2_analysis(DownSampledData)

# log2 fold change <0, TP 
ResultLTrue=Result[rownames(Result) %in% rownames(DAGs) & Result[,2]<0 & Result[,1]<0.05 ,]
summary(ResultLTrue[,2])

# log2 fold change <0, FP  
ResultLFalse=Result[!rownames(Result) %in% rownames(DAGs) & Result[,2]<0 & Result[,1]<0.05 ,]
summary(ResultLFalse[,2])

# log2 fold change >0, TP  
ResultBTrue=Result[rownames(Result) %in% rownames(DAGs) & Result[,2]>0 & Result[,1]<0.05 ,]
summary(ResultBTrue[,2])

# log2 fold change >0, FP  
ResultBFalse=Result[!rownames(Result) %in% rownames(DAGs) & Result[,2]>0 & Result[,1]<0.05 ,]
summary(ResultBFalse[,2])

#=================================================================================================================================

# Plot TP and FP in significance order 
PR=c()
for (i in 1:nrow(DownSampledData)) {
  
  if (BinTP[i]==0){
    PR[i]="FP"
  } else {
    PR[i]="TP"
  }
}
x=seq(1,nrow(DownSampledData), by=1)
y=c(rep(1, nrow(DownSampledData)))
HeatmapPR=data.frame(x,y, PR)

ggplot(HeatmapPR, aes(x=x, y=y, fill=PR)) +
  geom_tile(aes(fill = PR)) +
  #scale_fill_manual(values = c("red", "blue")) + 
  scale_fill_viridis_d(begin = 0, end = 1, alpha = 0.5) +  
  labs(title="test", x = "Group size", y = "Sequencing depth", fill = "AUC-values") 















# ============================================= Check for zeros ============================================= 
# Total number of zeros in the dataset
sum(DownSampledData==0)

# create vector containg the number of zeroes per gene
zerosPerGene=c()
for (i in 1:nrow(DownSampledData)) {
  zerosPerGene[i]=sum(DownSampledData[i,]==0)
}

summary(zerosPerGene)

# See if the gene with maximum zeros is among the introduced DAGs (if yes, returns at what index)
match(row.names(DownSampledData[match(92, zerosPerGene),]), rownames(DAGs))


