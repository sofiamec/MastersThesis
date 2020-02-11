# statistical analysis of DAGs in the resampled datasets

# Required "inputs"

# seed = selectedSeed 
# saveName = "Gut2"
# plotName = "Human Gut II"
# saveExpDesign = 
# plotExpDesign =


# load required packages
library(DESeq2)
library(edgeR)
library(DescTools)
library(ggplot2)

seed=selectedSeed  # set seed 

# resampled (and filtered) data without DAGS
# ResampData <- read.csv(file=sprintf("../../Intermediate/%s/%s/ResampData_seed%d.csv", saveName, saveExpDesign,seed), header = T, row.names = 1) 

# resampled data with DAGs
DagData <- read.csv(file=sprintf("../../Intermediate/%s/%s/downSampledData_seed%d.csv", saveName, saveExpDesign,seed), header = T, row.names = 1)
Dags <- read.csv(file=sprintf("../../Intermediate/%s/%s/DAGs_seed%d.csv", saveName, saveExpDesign, seed), header = T, row.names = 1)

#===================================================================================================================================
#=========================================== Functions ==============================================================================
#===================================================================================================================================

# Function for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

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

# Computing ROC-curves and AUC-values
# For the results from analysing DAGs in a dataset and the corresponding known DAGs,
# this function computes AUC-values and plots the ROC-curve.
# Inputs:   ResultsData = Results from DESeq2 or edgeR analysis. Ex: ResDESeq or ResEdge
#           Dags = the artificially introduced DAGs (known)
# Outputs:  ROC = a dataframe with the computed TPR- and FPR-values
#           AUCs = The computed AUC for the entire ROC-curve and for FPR-cutoff 0.05 and 0.10.
#           meanROC = the pieciwise mean
Compute_ROC_AUC = function(ResultsData, Dags, seed){
  
  TP<-rownames(Dags)
  nT=vector(mode = 'numeric' ,length = nrow(ResultsData)+1)
  nF=vector(mode = 'numeric' ,length = nrow(ResultsData)+1)
  for (i in 1:nrow(ResultsData)){
    if (match(rownames(ResultsData)[i],TP, nomatch = F)!=0){
      nT[i+1]=nT[i]+1
      nF[i+1]=nF[i]
    }
    else{
      nT[i+1]=nT[i]
      nF[i+1]=nF[i]+1
    }
  }
  
  nT<-nT[-c(1)]
  nF<-nF[-c(1)]
  
  TPR<-nT/nrow(Dags)
  FPR<-nF/(nrow(ResultsData)-nrow(Dags))
  
  # Beräkna på annat sätt!
  AUCtot <- AUC(FPR, TPR, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE), method = c("trapezoid"), absolutearea = FALSE, subdivisions = 100, na.rm = FALSE)
  AUC5 <- AUC(FPR, TPR, from = min(x, na.rm = TRUE), to = 0.05, method = c("trapezoid"), absolutearea = FALSE, subdivisions = 100, na.rm = FALSE)
  AUC10 <- AUC(FPR, TPR, from = min(x, na.rm = TRUE), to = 0.1, method = c("trapezoid"), absolutearea = FALSE, subdivisions = 100, na.rm = FALSE)
  
  # Ta ut TPR/FPR? vid fixerat värde
  
  ROCs <- data.frame(TPR,FPR,seed)
  ROCs[,2] <- round2(ROCs[,2], 3) # 3 is the number of decimals here
  
  meanROCs<-ddply(ROCs, "FPR", summarise,
                 N    = length(TPR),
                 mean = mean(TPR),
                 sd   = sd(TPR),
                 se   = sd / sqrt(N)
  )               
  AUCs <- data.frame(AUC5,AUC10,AUCtot,seed)
  return(list(ROCs, AUCs, meanROCs))
}

#===================================================================================================================================

######## DESeq2 ##########
ResDESeq=DESeq2_analysis(Data = DagData)
cat(sprintf("Number of significant genes with DESeq2 for %s: %d     (exp. design: %s)", plotName, sum(ResDESeq$padj<0.05),plotExpDesign))

# how many of the artificially introduced DAGs are among the significant genes
matchDESeq=c()
for (i in 1:nrow(Dags)) {
  matchDESeq[i]=sum(grepl(rownames(Dags)[i], rownames(ResDESeq[which(ResDESeq$padj<0.05),])))
}
cat(sprintf("Number of TP genes with DESeq2 for %s: %d              (exp. design: %s)", plotName, sum(matchDESeq), plotExpDesign)) 

rm(matchDESeq)

#===================================================================================================================================
# Computing ROC and AUC
# Plotting both deseq and edge (Lägg till detta i funktionen Compute_ROC_AUC när vi bestämt oss för edgeR eller DESeq!)
deseqROCAUC<-Compute_ROC_AUC(ResDESeq,Dags, seed)
ROCs <- data.frame(deseqROCAUC[[1]])
AUCs<- as.matrix(deseqROCAUC[[2]]) 
meanROCs<-as.matrix(deseqROCAUC[[3]])
rm(deseqROCAUC)