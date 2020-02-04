# statistical analysis of DAGs in the resampled datasets

# load required packages
library(DESeq2)
library(edgeR)

ResampData <- read.csv(file="../../Intermediate/ResampData.csv", header = T, row.names = 1) # read resampled dataset


#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================
# DESeq2-analysis
# This function uses DESeq2 to identfy DAGs in a dataset containing two groups
# Input: Data = the data to analyse
# Output: a dataframe containing the adjusted p-value and log2fold-change for each gene
DESeq2_analysis=function(Data){
  m=ncol(Data)/2                  # number of samples in each group in the dataset

  #DESeq2-analysis
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  ResultDESeq<-DESeq(CountsDataset)                                       # Perform analysis 
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(ResampData), Res$padj, Res$log2FoldChange)   # dataframe with genes, padj and log2fold-change
  colnames(Result) <- c("Gene", "padj", "log2foldchange")                 # name the columns
  
  return(Result)
}

#===================================================================================================================================


######## DESeq2 ##########
Result=DESeq2_analysis(ResampData)


######## edgeR ###########

edgeRUsersGuide(view = T)                    # user guide for edgeR
m=60                                         # number of samples

group <- factor(c(rep(1,m),rep(0,m)))        # grouping factor
design <- model.matrix(~group)               # design matrix
y <- DGEList(counts=ResampData, group=group) # combin dataset and grouping factor into a DGE-list
y <- estimateDisp(y, design, robust=TRUE)    # estimate the dispersion of the dataset
fit <- glmQLFit(y, design, robust = TRUE)    # fit the negative binomial GLM for each gene
qlf <- glmQLFTest(fit)                       # carries out the quasi-likelihood F-test

topTags(qlf)                                 # print the most significant DAGs
summary(decideTests(qlf))                    # summary of result





