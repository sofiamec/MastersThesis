# statistical analysis of DAGs in the resampled datasets

# load required packages
library(DESeq2)
library(edgeR)

seed=100  # set seed 

# resampled data without DAGS
ResampData <- read.csv(file=sprintf("../../Intermediate/ResampData_seed%d.csv",seed), header = T, row.names = 1) 

# resampled data with DAGs
DagData <- read.csv(file=sprintf("../../Intermediate/downSampledData_seed%d.csv",seed), header = T, row.names = 1)
Dags <- read.csv(file=sprintf("../../Intermediate/DAGs_seed%d.csv",seed), header = T, row.names = 1)

#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================
# DESeq2-analysis
# This function uses DESeq2 to identfy DAGs in a dataset containing two groups
# Input: Data = the data to analyse
# Output: a dataframe containing the adjusted p-value and log2fold-change for each gene, ordered with increading padj
DESeq2_analysis=function(Data){
  m=ncol(Data)/2                  # number of samples in each group in the dataset

  #DESeq2-analysis
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  ResultDESeq<-DESeq(CountsDataset)                                       # Perform analysis 
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(Data), Res$padj, Res$log2FoldChange)         # dataframe with genes, padj and log2fold-change
  Result <- data.frame(Result[,-1], row.names=Result[,1])                 # put first column (genes) as rowname
  colnames(Result) <- c("padj", "log2foldchange")                         # name the columns
  
  
  ResultSorted=Result[order(Result[,1]),]                                 # order with increasing padj
  
  return(ResultSorted)
}

#===================================================================================================================================
# edgeR-analysis
# This function uses edgeR to identfy DAGs in a dataset containing two groups
# Input: Data = the data to analyse
# Output: a dataframe containing the logFC, p-value and FDR for each gene, ordered with increasing FDR
edgeR_analysis=function(Data){
  m=ncol(Data)/2                                  # number of samples in each group in the dataset
  
  group <- factor(c(rep(1,m),rep(0,m)))           # grouping factor
  design <- model.matrix(~group)                  # design matrix
  y <- DGEList(counts=Data, group=group)          # combine dataset and grouping factor into a DGE-list
  y <- estimateDisp(y, design, robust=TRUE)       # estimate the dispersion of the dataset
  fit <- glmQLFit(y, design, robust = TRUE)       # fit the negative binomial GLM for each gene
  qlf <- glmQLFTest(fit)                          # carries out the quasi-likelihood F-test
  Out <- topTags(qlf, n = "Inf")$table[,c(1,4,5)] # print logFC, p-value and FDR

  OutSorted=Out[order(Out[,3]),]                  # order with increasing FDR
  
  return(OutSorted)
}
#===================================================================================================================================

################################## Comparison of methods ########################################
# Maybe remove when one method is chosen 

######## DESeq2 ##########
ResDESeq=DESeq2_analysis(Data = DagData)
sum(ResDESeq$padj<0.05)


# how many of the artificially introduced DAGs are among the significant genes
matchDESeq=c()
for (i in 1:nrow(Dags)) {
  matchDESeq[i]=sum(grepl(rownames(Dags)[i], rownames(ResDESeq[which(ResDESeq$padj<0.05),])))
}
sum(matchDESeq) 


######## edgeR ##########
ResEdge=edgeR_analysis(Data = DagData)
sum(ResEdge$FDR<0.05)

# how many of the artificially introduced DAGs are among the significant genes
matchEdge=c()
for (i in 1:nrow(Dags)) {
  matchEdge[i]=sum(grepl(rownames(Dags)[i], rownames(ResEdge[which(ResEdge$FDR<0.05),])))
}
sum(matchEdge)

