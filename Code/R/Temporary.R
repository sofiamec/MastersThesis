
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

#===================================================================================================================================
# Plot individual mean RoC-curves for certain experimental design
meanROCplot <- ggplot(data=meanROC2, aes(x=FPR, y=meanTPR)) +  theme_minimal() + 
  geom_ribbon(aes(ymin=(min), ymax=(max), fill="#22A88433"), alpha = 0.2) + 
  geom_line(aes(color="#22A88433")) + theme(legend.position = "none") +
  labs(title=sprintf("Mean ROC-curve for %s with effect %g", plotName,q), 
       subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
       x = "False Positive Rate", y = "True Positive Rate")+
  ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.6) +
  scale_colour_viridis_d(begin = 0.2, end = 0.6)

print(meanROCplot)

if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/IntermediatePlots/meanROC_%s.pdf", saveName, saveExpDesign)
  ggsave(filename = path_save, plot = meanROCplot, height = 5, width = 6)
  dev.off()
  print(meanROCplot)
  rm(path_save)
}
rm(meanROCplot)


