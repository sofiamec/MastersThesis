
### Exploration of Gut1
# check ordered sequencing depth for Gut1 
Gut1Depth=data.frame(colnames(Gut1), as.numeric(colSums(Gut1)))
gut1Order <- order(Gut1Depth[,2])
Gut1DepthSorted <- Gut1Depth[gut1Order,]  

gene_boxplot(Gut1GG,"Number of reads in Human Gut I", colorScale150, savePlot=F, saveName="Gut1")

sequencing_depth_histogram(Gut1GG, bins=30, "Sequencing Depth in Human Gut I", 
                           color1, "black", savePlot=F, saveName="Gut1")

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
######## edgeR ##########
ResEdge=edgeR_analysis(Data = DagData)
cat(sprintf("Number of significant genes with edgeR for %s: %d        (exp. design: %s)", plotName, sum(ResEdge$FDR<0.05)), plotExpDesign)

# how many of the artificially introduced DAGs are among the significant genes
matchEdge=c()
for (i in 1:nrow(Dags)) {
  matchEdge[i]=sum(grepl(rownames(Dags)[i], rownames(ResEdge[which(ResEdge$FDR<0.05),])))
}
cat(sprintf("Number of TP genes with edgeR for %s: %d                 (exp. design: %s)", plotName, sum(matchEdge)), plotExpDesign)


#===================================================================================================================================
# Plotting
edgeROCAUC<-Compute_ROC_AUC(ResEdge,Dags)
edgeROC <- data.frame(edgeROCAUC[[1]], "edgeR")
colnames(edgeROC)[4]<-"Dataset"




