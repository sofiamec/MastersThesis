#===================================================================================================================================
# ANALYSIS SECTION
#===================================================================================================================================

if (analysis=="DESeq"){
  
  ######## DESeq2 ##########
  ResDESeq=DESeq2_analysis(Data = DownSampledData)
  ResDAGsAnalysis<-ResDESeq
  rm(ResDESeq)
  
} else if (analysis=="OGLM") {
  
  ######## OGLM ##########
  ResOGLM=OGLM_analysis(Data = DownSampledData)
  ResDAGsAnalysis<-ResOGLM
  rm(ResOGLM)
  
} 

cat(sprintf("Significant genes for %s:    %d     (exp. design: %s)\n", saveName, sum(ResDAGsAnalysis[,2]<0.05, na.rm = T),plotExpDesign))

# how many of the artificially introduced DAGs are among the significant genes
matchDAGs=c()
for (i in 1:nrow(DAGs)) {
  matchDAGs[i]=sum(grepl(rownames(DAGs)[i], rownames(ResDAGsAnalysis[which(ResDAGsAnalysis[,2]<0.05),,drop=F])))
}
cat(sprintf("TP genes for %s out of %d:  %d     (exp. design: %s)\n\n", saveName, nrow(DAGs), sum(matchDAGs), plotExpDesign))

#===================================================================================================================================
# Computing ROC and AUC
ROCAUC<-Compute_ROC_AUC(ResDAGsAnalysis,DAGs, run, F)
ROCs <- data.frame(ROCAUC[[1]])
AUCs<- as.matrix(ROCAUC[[2]]) 
meanROCs<-as.matrix(ROCAUC[[3]])
genesFDRs<-as.matrix(ROCAUC[[4]])

if (runStrata==T){
  # Add corresponding strata to each gene in the analysed results
  ResStrata<-DataStrata[rownames(DataStrata) %in% rownames(ResDAGsAnalysis),]
  ResStrata<-ResStrata[rownames(ResDAGsAnalysis),]
  ResStrata<-data.frame(ResDAGsAnalysis, ResStrata, c(rep(0,nrow(ResDAGsAnalysis))))
  
  ROCsAbundance <- data.frame()
  AUCsAbundance <- data.frame()
  meanROCsAbundance <- data.frame()
  
  ROCsVariability <- data.frame()
  AUCsVariability <- data.frame()
  meanROCsVariability <- data.frame()
  
  for (k in 1:numberOfStrata) {
    DAGsStrataAbundance=DAGs[DAGs$AbundanceStrata==k,]
    DAGsStrataVariability=DAGs[DAGs$VariabilityStrata==k,]
    ResDAGsAnalysisAbundance=ResStrata[ResStrata$AbundanceStrata==k,]
    ResDAGsAnalysisVariability=ResStrata[ResStrata$VariabilityStrata==k,]
    
    deseqROCAUCAbundance<-Compute_ROC_AUC(ResDAGsAnalysisAbundance, DAGsStrataAbundance, run, T)
    deseqROCAUCVariability<-Compute_ROC_AUC(ResDAGsAnalysisVariability, DAGsStrataVariability, run, T)
    
    ROCsAbundance <- rbind(ROCsAbundance ,data.frame(deseqROCAUCAbundance[[1]], k))
    AUCsAbundance <- rbind(AUCsAbundance ,data.frame(deseqROCAUCAbundance[[2]],k))
    meanROCsAbundance <- rbind(meanROCsAbundance ,data.frame(deseqROCAUCAbundance[[3]],k))
    
    ROCsVariability <- rbind(ROCsVariability ,data.frame(deseqROCAUCVariability[[1]], k))
    AUCsVariability <- rbind(AUCsVariability ,data.frame(deseqROCAUCVariability[[2]],k))
    meanROCsVariability <- rbind(meanROCsVariability ,data.frame(deseqROCAUCVariability[[3]],k))
    
    rm(DAGsStrataAbundance, DAGsStrataVariability, ResDAGsAnalysisAbundance, ResDAGsAnalysisVariability,  deseqROCAUCAbundance, deseqROCAUCVariability)
  }
  
  colnames(AUCsAbundance)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","strata")
  colnames(AUCsVariability)<-c("AUC1", "AUC5" ,"AUCtot", "TPR1", "TPR5", "run","strata")
  colnames(meanROCsAbundance)<-c("FPR", "N", "meanTPR" ,"strata")
  colnames(meanROCsVariability)<-c("FPR", "N", "meanTPR" ,"strata")
  colnames(ROCsAbundance)<-c("TPR","FPR", "run","strata")
  colnames(ROCsVariability)<-c("TPR","FPR", "run","strata")
  
  rm(k, ResStrata)
}

# remove variables/datasets
rm(DownSampledData, matchDAGs, ROCAUC, ResDAGsAnalysis, DAGs, i)
