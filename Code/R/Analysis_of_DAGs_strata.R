# statistical analysis of DAGs in the resampled datasets

# Required "inputs"

# seed = selectedSeed 
# saveName = "Gut2"
# plotName = "Human Gut II"
# saveExpDesign = 
# plotExpDesign =
# run =

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
# Output: a dataframe containing the p-value and the adjusted p-value for each gene, ordered with increading p-values
DESeq2_analysis=function(Data){
  
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                     # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(Data), Res$pvalue, Res$padj)                   # dataframe with genes, their p-values and adjusted p-values
  ResultSorted=as.data.frame(Result[order(Result[,2]),])                    # order with increasing p-value
  ResultSorted <- data.frame(ResultSorted[,-1], row.names=ResultSorted[,1]) # put first column (genes) as rowname
  colnames(ResultSorted) <- c("p-value", "adjusted p-value")                # name the columns
  
  return(ResultSorted)
}

#===================================================================================================================================
# Computing ROC-curves and AUC-values
# For the results from analysing DAGs in a dataset and the corresponding known DAGs,
# this function computes AUC-values and plots the ROC-curve.
# Inputs:   ResultsData = Results from DESeq2 or edgeR analysis. Ex: ResDESeq or ResEdge
#           DAGs = the artificially introduced DAGs (known)
# Outputs:  ROC = a dataframe with the computed TPR- and FPR-values
#           AUCs = The computed AUC for the entire ROC-curve and for FPR-cutoff 0.05 and 0.10.
#           meanROC = the pieciwise mean
Compute_ROC_AUC = function(ResultsData, trueTP, run){
  
  TP<-rownames(trueTP)
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
  
  TPR<-nT/nrow(trueTP)
  FPR<-nF/(nrow(ResultsData)-nrow(trueTP))
  
  # Compute AUC and TPR at certain FPR
  AUC5<-trapz(FPR[FPR<=0.05],TPR[FPR<=0.05])/max(FPR[FPR<=0.05])
  AUC1<-trapz(FPR[FPR<=0.01],TPR[FPR<=0.01])/max(FPR[FPR<=0.01])
  AUCtot<-trapz(FPR,TPR)
  TPR5<-max(TPR[FPR<=0.05])
  TPR1<-max(TPR[FPR<=0.01])
  AUCs <- data.frame(AUC1,AUC5,AUCtot,TPR1,TPR5,run)
  
  rm(AUC1,AUC5,AUCtot, TPR1,TPR5)
  
  ROCs <- data.frame(TPR,FPR,run)
  ROCs2 <- ROCs
  ROCs2[,2] <- round2(ROCs2[,2], 3) # 3 is the number of decimals here
  
  meanROCs<-ddply(ROCs2, "FPR", summarise,
                 N    = length(TPR),
                 mean = mean(TPR))
  
  return(list(ROCs, AUCs, meanROCs))
}

#===================================================================================================================================

######## DESeq2 ##########
ResDESeq=DESeq2_analysis(Data = DownSampledData)
cat(sprintf("Significant genes for %s:    %d     (exp. design: %s)\n", saveName, sum(ResDESeq[,2]<0.05, na.rm = T),plotExpDesign))

# how many of the artificially introduced DAGs are among the significant genes
matchDESeq=c()
for (i in 1:nrow(DAGs)) {
  matchDESeq[i]=sum(grepl(rownames(DAGs)[i], rownames(ResDESeq[which(ResDESeq[,2]<0.05),,drop=F])))
}
cat(sprintf("TP genes for %s out of %d:  %d     (exp. design: %s)\n\n", saveName, nrow(DAGs), sum(matchDESeq), plotExpDesign))

#===================================================================================================================================
# Computing ROC and AUC
# Plotting both deseq and edge (Lägg till detta i funktionen Compute_ROC_AUC när vi bestämt oss för edgeR eller DESeq!)
deseqROCAUC<-Compute_ROC_AUC(ResDESeq,DAGs, run)
ROCs <- data.frame(deseqROCAUC[[1]])
AUCs<- as.matrix(deseqROCAUC[[2]]) 
meanROCs<-as.matrix(deseqROCAUC[[3]])

# Add corresponding strata to each gene in the analysed results
ResStrata<-DataStrata[rownames(DataStrata) %in% rownames(ResDESeq),]
ResStrata<-ResStrata[rownames(ResDESeq),]
ResStrata<-data.frame(ResDESeq, ResStrata, c(rep(0,nrow(ResDESeq))))
#ResStrata[rownames(ResStrata) %in% rownames(DAGs),7]=1 #Adding notation for which are TP



for (k in 1:length(numberOfStrata)) {
  DAGsStrataAbundance=DAGs[DAGs$AbundanceStrata==k,]
  DAGsStrataVariability=DAGs[DAGs$VariabilityStrata==k,]
  ResDESeqAbundance=ResStrata[ResStrata$AbundanceStrata==k,]
  ResDESeqVariability=ResStrata[ResStrata$VariabilityStrata==k,]
  
  deseqROCAUCAbundance<-Compute_ROC_AUC(ResDESeqAbundance, DAGsStrataAbundance, run)
  deseqROCAUCVariability<-Compute_ROC_AUC(ResDESeqVariability, DAGsStrataVariability, run)
  
  ROCsAbundance <- data.frame(deseqROCAUCAbundance[[1]], k)
  AUCsAbundance <- data.frame(deseqROCAUCAbundance[[2]],k) 
  meanROCsAbundance <- data.frame(deseqROCAUCAbundance[[3]],k)
  
  ROCsVariability <- data.frame(deseqROCAUCVariability[[1]], k)
  AUCsVariability <- data.frame(deseqROCAUCVariability[[2]],k) 
  meanROCsVariability <- data.frame(deseqROCAUCVariability[[3]],k)  
}

colnames(AUCsAbundance)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","Strata")
colnames(AUCsVariability)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","Strata")
colnames(meanROCsAbundance)<-c("FPR", "N", "meanTPR" ,"Strata")
colnames(meanROCsVariability)<-c("FPR", "N", "meanTPR" ,"Strata")
colnames(ROCsAbundance)<-c("TPR","FPR", "run","Strata")
colnames(ROCsVariability)<-c("TPR","FPR", "run","Strata")



# remove variables/datasets
rm(DownSampledData, matchDESeq, deseqROCAUC, ResDESeq, DAGs, i)
