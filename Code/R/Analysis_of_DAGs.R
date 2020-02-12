# statistical analysis of DAGs in the resampled datasets

# load required packages
library(DESeq2)
library(ggplot2)
library(pracma)

#===================================================================================================================================
## Selecting parameters and data:

saveName = "Gut2" # Choose dataset. Ex: "Gut2" or "Marine"
m = 60        # Number of samples in each group (total nr samples = 2*m)
d = 10000    # Desired sequencing depth per sample. It will not be exct
q = 2         # Fold-change for downsampling
f = 0.10      # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
seed=1  
savePlot=F

{ # Quickly gives the case the correct names
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
  saveExpDesign = sprintf("m%d_d%d%s_q%d_f%d", m, dD, prefix, q, f*100)
  plotExpDesign = sprintf("m=%d, d=%d%s, q=%d, f=%d%%",m,dD,prefix,q,f*100)

rm(dD,prefix)

# Create folder for certain case if it doesn't exist
if (!dir.exists(sprintf("../../Result/%s/%s", saveName, saveExpDesign))){
  dir.create(file.path("../../Result", sprintf("%s", saveName), sprintf("%s", saveExpDesign)), recursive = T)
}
  
# resampled data with DAGs
DagData <- read.csv(file=sprintf("../../Intermediate/%s/%s/DownSampledData_seed%d.csv",saveName,saveExpDesign,seed), header = T, row.names = 1)
Dags <- read.csv(file=sprintf("../../Intermediate/%s/%s/DAGs_seed%d.csv", saveName, saveExpDesign,seed), header = T, row.names = 1)
}
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
# Output: a dataframe containing the p-value and log2fold-change for each gene, ordered with increading p-values
DESeq2_analysis=function(Data){
  m=ncol(Data)/2                  # number of samples in each group in the dataset

  #DESeq2-analysis
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                     # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(Data), Res$pvalue, Res$log2FoldChange)       # dataframe with genes, p-values and log2fold-change
  Result <- data.frame(Result[,-1], row.names=Result[,1])                 # put first column (genes) as rowname
  colnames(Result) <- c("p-value", "log2foldchange")                      # name the columns
  
  
  ResultSorted=Result[order(Result[,1]),]                                 # order with increasing p-value
  
  return(ResultSorted)
}

#===================================================================================================================================

# Computing ROC-curves and AUC-values
# For the results from analysing DAGs in a dataset and the corresponding known DAGs,
# this function computes AUC-values and plots the ROC-curve.
# Inputs:   ResultsData = Results from DESeq2 or edgeR analysis. Ex: ResDESeq or ResEdge
#           Dags = the artificially introduced DAGs (known)
#           seed = the number of the selected seed
#           savePlot = TRUE/FALSE indicates wether or not to save the ROC-plot
# Outputs:  ROC = a dataframe with the computed TPR- and FPR-values
#           AUCs = The computed AUC for the entire ROC-curve and for FPR-cutoff 0.05 and 0.10.
#           meanROC = the pieciwise mean
Compute_ROC_AUC = function(ResultsData, Dags, seed, plotExpDesign, plotName, saveName, saveExpDesign, savePlot){
  
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
  
  # Compute AUC and TPR at certain FPR
  AUC5<-trapz(FPR[FPR<=0.05],TPR[FPR<=0.05])/max(FPR[FPR<=0.05])
  AUC10<-trapz(FPR[FPR<=0.1],TPR[FPR<=0.1])/max(FPR[FPR<=0.1])
  AUCtot<-trapz(FPR,TPR)
  TPR5<-max(TPR[FPR<=0.05])
  TPR10<-max(TPR[FPR<=0.1])
  AUCs <- data.frame(AUC5,AUC10,AUCtot,TPR5,TPR10,seed)
  
  rm(AUC5,AUC10,AUCtot, TPR5,TPR10)
  
  ROCs <- data.frame(TPR,FPR,seed)
  ROCs2<-ROCs
  ROC2s[,2] <- round2(ROCs2[,2], 3) # 3 is the number of decimals here
  
  meanROCs<-ddply(ROCs2, "FPR", summarise,
                  N    = length(TPR),
                  mean = mean(TPR),
                  sd   = sd(TPR),
                  se   = sd / sqrt(N) )
  rm(ROCs2)
  
  ROCplot <- ggplot(data=ROCs, aes(x=FPR, y=TPR)) +  geom_line() + 
    theme(plot.title = element_text(hjust = 0.5)) +  theme_minimal() + 
    scale_color_manual(values=c('#7FCDBB','#225EA8')) +
    labs(title=sprintf("ROC-curves for analysis of %s seed %d", plotName, seed), 
         subtitle = sprintf("Experimental design: %s    (Seed %d)", plotExpDesign, seed),
         x = "False Positive Rate", y = "True Positive Rate")
  
  print(ROCplot)
  
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/%s/ROC_seed%d.pdf",saveName, saveExpDesign, seed)
    ggsave(filename = path_save2, plot = p2, height = 5, width = 6)
    dev.off()
    print(p2)
  }
  return(list(ROCs, AUCs, meanROCs))
}

#===================================================================================================================================

######## DESeq2 ##########
ResDESeq=DESeq2_analysis(Data = DagData)
cat(sprintf("Number of significant genes with DESeq2 for %s: %d     (exp. design: %s)\n", plotName, sum(ResDESeq$padj<0.05),plotExpDesign))

# how many of the artificially introduced DAGs are among the significant genes
matchDESeq=c()
for (i in 1:nrow(Dags)) {
  matchDESeq[i]=sum(grepl(rownames(Dags)[i], rownames(ResDESeq[which(ResDESeq$padj<0.05),])))
}
cat(sprintf("Number of TP genes with DESeq2 for %s: %d              (exp. design: %s)\n", plotName, sum(matchDESeq), plotExpDesign))

#===================================================================================================================================
# Computing ROC and AUC
deseqROCAUC<-Compute_ROC_AUC(ResDESeq,Dags, seed, plotExpDesign, plotName, saveName, SaveExpDesign, savePlot)
ROCs <- data.frame(deseqROCAUC[[1]])
AUCs<- as.matrix(deseqROCAUC[[2]]) 
meanROCs<-as.matrix(deseqROCAUC[[3]])
