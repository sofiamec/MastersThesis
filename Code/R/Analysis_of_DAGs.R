#===================================================================================================================================
# ANALYSIS SECTION
#===================================================================================================================================
# statistical analysis of DAGs in the resampled datasets

# Required "inputs"

# seed = selectedSeed 
# saveName = "Gut2"
# plotName = "Human Gut II"
# saveExpDesign = 
# plotExpDesign =
# run =

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
deseqROCAUC<-Compute_ROC_AUC(ResDESeq,DAGs, run, F)
ROCs <- data.frame(deseqROCAUC[[1]])
AUCs<- as.matrix(deseqROCAUC[[2]]) 
meanROCs<-as.matrix(deseqROCAUC[[3]])
genesFDRs<-as.matrix(deseqROCAUC[[4]])

if (runStrata==T){
  # Add corresponding strata to each gene in the analysed results
  ResStrata<-DataStrata[rownames(DataStrata) %in% rownames(ResDESeq),]
  ResStrata<-ResStrata[rownames(ResDESeq),]
  ResStrata<-data.frame(ResDESeq, ResStrata, c(rep(0,nrow(ResDESeq))))
  
  ROCsAbundance <- data.frame()
  #AUCsAbundance <- data.frame()
  meanROCsAbundance <- data.frame()
  
  ROCsVariability <- data.frame()
  #AUCsVariability <- data.frame()
  meanROCsVariability <- data.frame()
  
  for (k in 1:numberOfStrata) {
    DAGsStrataAbundance=DAGs[DAGs$AbundanceStrata==k,]
    DAGsStrataVariability=DAGs[DAGs$VariabilityStrata==k,]
    ResDESeqAbundance=ResStrata[ResStrata$AbundanceStrata==k,]
    ResDESeqVariability=ResStrata[ResStrata$VariabilityStrata==k,]
    
    deseqROCAUCAbundance<-Compute_ROC_AUC(ResDESeqAbundance, DAGsStrataAbundance, run, T)
    deseqROCAUCVariability<-Compute_ROC_AUC(ResDESeqVariability, DAGsStrataVariability, run, T)
    
    ROCsAbundance <- rbind(ROCsAbundance ,data.frame(deseqROCAUCAbundance[[1]], k))
    #AUCsAbundance <- rbind(AUCsAbundance ,data.frame(deseqROCAUCAbundance[[2]],k))
    meanROCsAbundance <- rbind(meanROCsAbundance ,data.frame(deseqROCAUCAbundance[[3]],k))
    
    ROCsVariability <- rbind(ROCsVariability ,data.frame(deseqROCAUCVariability[[1]], k))
    #AUCsVariability <- rbind(AUCsVariability ,data.frame(deseqROCAUCVariability[[2]],k))
    meanROCsVariability <- rbind(meanROCsVariability ,data.frame(deseqROCAUCVariability[[3]],k))
    
    rm(DAGsStrataAbundance, DAGsStrataVariability, ResDESeqAbundance, ResDESeqVariability,  deseqROCAUCAbundance, deseqROCAUCVariability)
  }
  
  #colnames(AUCsAbundance)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","strata")
  #colnames(AUCsVariability)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","strata")
  colnames(meanROCsAbundance)<-c("FPR", "N", "meanTPR" ,"strata")
  colnames(meanROCsVariability)<-c("FPR", "N", "meanTPR" ,"strata")
  colnames(ROCsAbundance)<-c("TPR","FPR", "run","strata")
  colnames(ROCsVariability)<-c("TPR","FPR", "run","strata")
  
  rm(k, ResStrata)
}

# remove variables/datasets
rm(DownSampledData, matchDESeq, deseqROCAUC, ResDESeq, DAGs, i)
