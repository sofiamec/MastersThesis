#===================================================================================================================================
# FUNCTION AND LIBRARIES SECTION
#===================================================================================================================================
# The listed libraries and functions in this script are required for running Run_entire_analysis

#===========================================================================================================================================
#=========================================== Global Libraries ==============================================================================
## Loading Libraries:

library(plyr)
library(viridis)
library(DESeq2)
library(ggplot2)
library(pracma)
library(readxl)
library(RColorBrewer)

#===========================================================================================================================================
#=========================================== Global Functions ==============================================================================

## FUNCTION for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
    z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
    return(z)
}

#===========================================================================================================================================
#=========================================== Setup Functions ===============================================================================

## FUNCTION to Remove low counts 
# For a given dataset, this function removes genes with low counts (>75 % or an average count <3).
# input: Data = the data to remove genes from
# output: 
remove_low_counts=function(Data){
  
  if (saveName!="Resistance"){
    a=rowMeans(Data)<3
    b=vector()
    r=vector()
    for (i in 1:nrow(Data)) {
      b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
      r[i]<-a[i]+b[i]
    }
    FilteredData=Data[r==0,]
    
  } else if (saveName=="Resistance") {
    # ?ndra eventuellt senare!
    a=rowMeans(Data)<6
    FilteredData=Data[a==F,]
  }
  return(FilteredData)
}

## FUNCTION for creating strata
# DESeq2-analysis for original data 
# This function uses DESeq2 to estimate the mean count and dispersion for genes in a dataset 
# (after normalising them based on sequencing depth). This is then used for creating strata.
#
# Input:  Data = the data to analyse (Gut2 or Marine)
#         numberOfStrata = the number of strata/levels used when dividing the data
# Output: Result = a dataframe containing the mean values and the estimated dispersion for each gene 
#                  as well as the corresponding abundance- and variability-strata for each gene
DESeq2_for_strata=function(Data,numberOfStrata){
  
  DesignMatrix <- data.frame(group=factor(rep(1,(ncol(Data)))))                        # all samples belong to the same group 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data,DesignMatrix, design=~1) # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                           # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq,independentFiltering=FALSE,cooksCutoff=FALSE)         # extract results
  
  Result <- data.frame(rownames(Data), Res$baseMean, dispersions(ResultDESeq))  # dataframe with genes, their mean values (after normalization) and dispersion estimates
  Result <- data.frame(Result[,-1], row.names=Result[,1])                       # put first column (genes) as rowname
  
  # create a vector for assigning smallest to highest strata
  n=numberOfStrata
  repVector<-as.integer(cumsum(c(rep(nrow(Data)/n,n))))
  repVector<-repVector-c(0,repVector)[-(length(repVector)+1)]
  strataVector<-rep(1:n,repVector)
  
  # sort by abundance and assign abundance-strata
  Result<-Result[order(Result[,1]),] 
  Result<-cbind(Result,strataVector)
  
  # sort by variability and assign variability-strata
  Result<-Result[order(Result[,2]),]
  Result<-cbind(Result,strataVector)
  
  # rename columns, factorise strata and sort genes by original order
  colnames(Result)<-c("BaseMean", "Estimated dispersion", "AbundanceStrata", "VariabilityStrata")
  Result$AbundanceStrata<-as.factor(Result$AbundanceStrata)
  Result$VariabilityStrata<-as.factor(Result$VariabilityStrata)
  Result<-Result[rownames(Data),]
  
  rm(n, numberOfStrata, repVector,strataVector)
  return(Result)
}


#===========================================================================================================================================
#=========================================== Resample Functions ============================================================================

# FUNCTION for Resampling
# This function resamples a new datasets with m*2 samples and d reads in each sample (depth) 
# input arguments:
# Data = the data to sample from
# m = number of samples in the new datasets
# d = sequencing depth for each sample in the new datasets
# output:  the resampled data in a large dataframe containing m*2 groups
resample = function(Data, m, d){
  sampleVector=sample(ncol(Data), 2*m)                      # vector w. the columnnumber of the sampled samples for both datasets
  DataNew=data.frame(row.names(Data), stringsAsFactors = F) # dataframe to put the resampled data in
  
  for (i in 1:(2*m)){
    readList <- rep(row.names(Data), times=Data[,sampleVector[i]])                # vector containing each read as one entry 
    sampledReads <- sample(readList, size=d)                                      # vector with d resampled reads from "readList"
    sampledVector <- as.data.frame(table(sampledReads), stringsAsFactors = FALSE) # assemble vector "sampledReads" into dataframe
    colnames(sampledVector) <- c("Gene",colnames(Data[sampleVector[i]]))          # name the columns according to the sample-name  
    DataNew <- merge(DataNew, sampledVector, by.x = 1, by.y = 1, all.x = T)       # insert "sampledVector" to "DataNew"
  }
  
  DataNew <- data.frame(row.names=DataNew[,1], DataNew[,-1]) # put first column (genes) as rownames  
  DataNew[is.na(DataNew)] <- as.integer(0)                   # set all "NA" to 0 (as integers, since DESeq2 require integer counts)
  
  return(DataNew) 
}

# FUNCTION for Summary of low counts
# For a given dataset, this function computes genes which should be removed when filtering.
# It also prints how many genes will be removed respectively kept.
# input: Data = the data to be filtered
# output: r = a logical vector with 1 if a gene should be removed, and 0 if it should be kept. 
compute_low_counts=function(Data){
  
  a=rowMeans(Data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(Data)) {
    b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
    r[i]<-a[i]+b[i]
  }
  cat(sprintf("Genes with low counts:         %s\n", sum(r!=0)))
  cat(sprintf("Genes with acceptable counts:  %s\n", sum(r==0)))
  return(r)
}

# FUNCTION for introducing DAGs
# For a resampled dataset (including both groups), this function introduces DAGs by
# downsampling a given fraction of genes. The DAGs will be balanced in the two groups.
# Input: Data = the resampled dataset including both groups (dataset 1 and 2)
#        q = the fold change. Ex. 10. Will result in relative abundance between datasets.
#        f = the desired (total) fraction of genes to be downsampled (the output will not follow this exactly)
# Outut: downSampledData = the new dataset with introduced DAGs. Analysis should be performed on this dataset
#        DAGs = a matrix with an overview of which genes have been downsampled in which dataset with the given q
introducing_DAGs = function(Data, q, f){
  
  downSampledData = Data
  if (saveName=="Resistance"){
    downSampledData = downSampledData[!rownames(downSampledData) %in% "Non-resistance",]
  }
  nDAGs = 2 * round(f*nrow(downSampledData)/2) # the total number of genes to be downsampled if we don't allow unbalanced DAGs
  #nDAGs = trunc(f*nrow(downSampledData)+0.5) # the total number of genes to be downsampled if we allow unbalanced DAGs
  
  # Creating empty matrix for overview of DAGs
  DAGs = matrix(ncol = 2, nrow = nrow(downSampledData)) 
  colnames(DAGs) = c(sprintf("Sample 1 to %d", ncol(Data)/2),sprintf("Sample %d to %d", ncol(Data)/2+1, ncol(Data)))
  rownames(DAGs) <- rownames(downSampledData)
  # Selecting random genes
  randomGenes <- sample(nrow(downSampledData),nDAGs) # Selects n random genes in the dataset which will be downsampled. 
  if (nDAGs==1|| nDAGs==0) {
    cat(print("Too few DAGs introduced"))
    nDAGs=2
  }
  # if we don't allow unbalanced DAGs
  rG1 <- randomGenes[1:(nDAGs/2)] # will be downsampled in dataset 1 
  rG2 <- randomGenes[(nDAGs/2+1):nDAGs] # will be downsampled in dataset 2
  
  for (gene in rG1) {
    for (sample in 1:(ncol(downSampledData)/2)) {
      downSampledData[gene,sample] <- rbinom(n = 1 ,size = downSampledData[gene,sample] ,prob = 1/q)
    }
    DAGs[gene, 1] <- -q
  }
  
  for (gene in rG2) {
    for (sample in (ncol(downSampledData)/2+1):(ncol(downSampledData))) {
      downSampledData[gene,sample] <- rbinom(n = 1 ,size = downSampledData[gene,sample] ,prob = 1/q) 
      DAGs[gene, 2] <- -q
    }
  }
  
  if (saveName=="Resistance"){
    downSampledData = rbind(Data[rownames(Data) %in% "Non-resistance",], downSampledData)
    downSampledData <- downSampledData[order(rownames(downSampledData)),]
  }
  
  DAGs<-DAGs[rowSums(DAGs, na.rm=T)!=0,]
  return(list(downSampledData,DAGs))

}

#===========================================================================================================================================
#=========================================== Analysis of DAGs Functions ====================================================================

## FUNCTION for DESeq2-analysis
# This function uses DESeq2 to identfy DAGs in a dataset containing two groups
# Input: Data = the data to analyse
# Output: a dataframe containing the p-value and the adjusted p-value for each gene, ordered with increading p-values
DESeq2_analysis=function(Data){
  
  DesignMatrix <- data.frame(group=factor(c(rep(1,m),rep(0,m))))          # define the different groups 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data, 
                                        DesignMatrix, design=~group)      # combine design matrix and data into a dataset
  
  if (saveName == "Resistance"){                                          # create special case for Resisance since it contains few genes 
    if(m==3 | m==5 | d==10000){                                                      # all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value, and so the standard curve fitting techniques will not work.
                                                                          # Perform DESeqs steps separatley:
      CountsDataset <- estimateSizeFactors(CountsDataset) 
      CountsDataset <- estimateDispersionsGeneEst(CountsDataset)
      dispersions(CountsDataset) <- mcols(CountsDataset)$dispGeneEst      # Use the gene-wise estimates as final estimates
      ResultDESeq <- nbinomWaldTest(CountsDataset)
    } else {
      ResultDESeq<-suppressMessages(DESeq(CountsDataset, fitType = "mean")) # Perform analysis with non-default fit-Type (suppress messages from it)  
    }
  } else {
    ResultDESeq<-suppressMessages(DESeq(CountsDataset))                     # Perform analysis (suppress messages from it)  
  }
  Res=results(ResultDESeq, independentFiltering=FALSE, cooksCutoff=FALSE) # extract results
  
  Result=data.frame(rownames(Data), Res$pvalue, Res$padj)                   # dataframe with genes, their p-values and adjusted p-values
  ResultSorted=as.data.frame(Result[order(Result[,2]),])                    # order with increasing p-value
  ResultSorted <- data.frame(ResultSorted[,-1], row.names=ResultSorted[,1]) # put first column (genes) as rowname
  colnames(ResultSorted) <- c("p-value", "adjusted p-value")                # name the columns
  
  return(ResultSorted)
}

## FUNCTION for OGLM-analysis
# This function uses OGLM and ANOVA to identfy DAGs in a dataset containing two groups
# Input:  Data = the data to analyse (DownsampledData)
# Output: ResOGLM = a dataframe containing the p-value and the adjusted p-value for each gene, ordered with increading p-values
OGLM_analysis<-function(Data){
  
  # Transposing and adding Sample group 
  ModData<-data.frame(c(rep(1,m),rep(0,m)),t(Data)) # DO NOT OPEN/VIEW TestA !!
  colnames(ModData)[1]<-"SampleGroup"
  
  logTotSeq<-log(colSums(Data))
  GeneNames<-c("Sample group",rownames(Data))
  Results <- data.frame(Gene=numeric(0), pValue=numeric(0))
  
  # Loopin over all genes and comparing models with and without group-covariate
  for (gene in 2:ncol(ModData)) {
    
    Model_Group<-glm(ModData[,gene] ~ SampleGroup + offset(logTotSeq), family = quasipoisson(link = "log"), data = ModData)
    Model_NoGroup<-glm(ModData[,gene] ~ offset(logTotSeq), family = quasipoisson(link = "log"), data = ModData)
    
    ANOVAres<-anova(Model_Group, Model_NoGroup, test = "F")
    Results<-rbind(Results,data.frame(GeneNames[gene],ANOVAres$`Pr(>F)`[[2]]))
  }
  
  Results<-data.frame(row.names = Results[,1], Results[,2],p.adjust(Results[,2], method = "BH"))
  colnames(Results)<-c("p-value", "adjusted p-value")
  ResOGLM<-Results[order(Results[,1]),]
  
  return(ResOGLM)
}


## FUNCTION for t-test-analysis
# This function uses t-test to identfy DAGs in a dataset containing two groups
# Input:  Data = the data to analyse (DownsampledData)
# Output: ResTTest = a dataframe containing the p-value and the adjusted p-value for each gene, ordered with increading p-values

ttest_analysis=function(Data){
  GeneNames<-rownames(Data)                                          
  Results <- data.frame(Gene=numeric(0), pValue=numeric(0))     # Dataframe to put results in
  
  for (i in 1:nrow(Data)) {                                     # loop over all genes, make one t-test for each gene
    
    # If both groups have variance=0 the p-value is set to NA
    if(var(as.numeric(sqrt(Data[i,1:m])))==0 & var(as.numeric(sqrt(Data[i,(m+1):(2*m)])))==0){
      Results=rbind(Results, data.frame(Gene=GeneNames[i],pValue=NA))
    }
    # If not, calculate p-values with t-test
    else {
      Res=t.test(sqrt(Data[i,1:m]),sqrt(Data[i,(m+1):(2*m)]))
      Results=rbind(Results, data.frame(Gene=GeneNames[i],pValue=Res$p.value)) 
    }
  }
  
  Results<-data.frame(row.names = Results[,1], Results[,2], 
                      p.adjust(Results[,2], method = "BH"))       # calculate adjusted p-values
  colnames(Results)<-c("p-value", "adjusted p-value")             # structure results
  ResTTest<-Results[order(Results[,1]),]                          # order results by increasing p-values
  
  return(ResTTest)
}

## FUNCTION for computing ROC-curves and AUC-values
# For the results from analysing DAGs in a dataset and the corresponding known DAGs,
# this function computes AUC-values and plots the ROC-curve.
# Inputs:   ResultsData = Results from DESeq2 or edgeR analysis. Ex: ResDESeq or ResEdge
#           DAGs = the artificially introduced DAGs (known)
# Outputs:  ROC = a dataframe with the computed TPR- and FPR-values
#           AUCs = The computed AUC for the entire ROC-curve and for FPR-cutoff 0.05 and 0.10.
#           meanROC = the pieciwise mean
Compute_ROC_AUC = function(ResultsData, trueTP, run, computeStrata){
  
  if (saveName=="Resistance"){
    ResultsData<-ResultsData[!rownames(ResultsData) %in% "Non-resistance",]
  }
  
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
  
  # with this muted (0,0) is included
  #nT<-nT[-c(1)]
  #nF<-nF[-c(1)]
  
  TPR<-nT/nrow(trueTP)
  FPR<-nF/(nrow(ResultsData)-nrow(trueTP))
  
  # Compute AUC and TPR at certain FPR
  AUC5<-trapz(FPR[FPR<=0.05],TPR[FPR<=0.05])/max(FPR[FPR<=0.05])
  AUC1<-trapz(FPR[FPR<=0.01],TPR[FPR<=0.01])/max(FPR[FPR<=0.01])
  AUCtot<-trapz(FPR,TPR)
  TPR5<-max(TPR[FPR<=0.05])
  TPR1<-max(TPR[FPR<=0.01])
  AUCs <- data.frame(AUC1,AUC5,AUCtot,TPR1,TPR5,run)
  
  # Compute nrTP nrFP at FDR 0.05
  cutoffID<-sum(ResultsData[,2]<0.05, na.rm = T)
  if (cutoffID==0){
    genesFDR<-data.frame(0,0,NA)
  } else {
    genesFDR<-data.frame(nT[cutoffID+1], nF[cutoffID+1],(nF[cutoffID+1]/(nF[cutoffID+1]+nT[cutoffID+1]))) 
  }
  colnames(genesFDR)<-c("NumberOfTP","NumberOfFP","TrueFDR")
  
  rm(AUC1,AUC5,AUCtot, TPR1,TPR5)
  
  ROCs <- data.frame(TPR,FPR,run)
  ROCs2 <- ROCs
  
  if (computeStrata==F){
    ROCs2[ROCs2[,2]<=0.05,2] <- round2(ROCs2[ROCs2[,2]<=0.05,2], 5) # 5 is the number of decimals here
    ROCs2[ROCs2[,2]>0.05,2] <- round2(ROCs2[ROCs2[,2]>0.05,2], 3) # 3 is the number of decimals here
  } else if (computeStrata==T){
    ROCs2[,2] <- round2(ROCs2[,2], 2) # 2 is the number of decimals here
  }
  
  meanROCs<-ddply(ROCs2, "FPR", summarise,
                  N    = length(TPR),
                  mean = mean(TPR))
  
  
  return(list(ROCs, AUCs, meanROCs, genesFDR))
}

#===========================================================================================================================================
#=========================================== Run Entire Functions ==========================================================================

## FUNCTION for plotting ROC-curves for all runs of an experimental design
# Inputs: ROCData = a combined matrix for all runs with exact FPR- and TPR-values as well as corresponding run
#         the function also needs other parameters, but it works as long as it is called on within the loop over d (after run-loop)
# Output  a plot is always printed, and if savePlot=T it will be saved
# Plot individual ROC-plots
individual_ROC_plot <- function(ROCData){
  ROCplot <- ggplot(data=ROCData, aes(x=FPR, y=TPR, group=run)) +  geom_line(aes(color=run)) + 
    scale_color_viridis(begin = 0, end = 0.85, discrete=TRUE) +
    labs(title=sprintf("ROC curves for %s with effect %g", plotName, q), 
         subtitle = sprintf("Experimental design: %s", plotExpDesign),
         colour="run", x = "False Positive Rate", y = "True Positive Rate") +
    ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))
  print(ROCplot)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/IntermediatePlots/individualROCs_%s.pdf", saveName, saveExpDesign)
    ggsave(filename = path_save, plot = ROCplot, height = 5, width = 6)
    dev.off()
    print(ROCplot)
    rm(path_save)
  }
  rm(ROCplot)
}

## FUNCTION for plotting AUC- and TPR-heatmaps
# Inputs:   variable = HeatmapData$AUCtot,HeatmapData$TPR1 etc
#           variableName = "total AUC-values", "TPR-values at FPR 0.01" etc
#           fillName = "AUC-values" or "TPR-values"
#           variableSave = "AUCtot", "TPR1" etc
# Outputs:  a plotted heatmap which is saved if savePlot==T
plot_heatmaps<-function(variable,variableName, fillName, variableSave){

  if (variableSave!="FDR"){
    fillCondition=variable 
    fillScale=scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) 
    
  } else if (variableSave=="FDR"){
    variableText = round2(variable, 2)
    cutOffs <-rev(c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
    for (i in 1:length(cutOffs)){
      variableText[variable <= cutOffs[i]] = sprintf("\u2264 %g", cutOffs[i])
    }
    fillCondition = variableText
    fillScale = scale_fill_manual(values = alpha(c("\u2264 1" = "#9E0142", "\u2264 0.9" = "#8F1257", "\u2264 0.8" = "#81236C", "\u2264 0.7" = "#733582", "\u2264 0.6" = "#654697", "\u2264 0.5" = "#5955A4", "\u2264 0.4" = "#4F62AB", "\u2264 0.3" = "#456EB1", "\u2264 0.2" = "#3B7BB7", "\u2264 0.1" = "#3288BD", "\u2264 0.05" = "#00BC77"), .6), na.value="grey50")
  }

  heatmap <-ggplot(HeatmapData, aes(x=m, y=d, fill=fillCondition)) +
    geom_tile(aes(fill = fillCondition)) + geom_text(aes(label = round2(variable, 2), fontface=md)) +
    fillScale + scale_y_discrete(labels = c("10000" = "10 k", "1e+05" = "100 k", "5e+05" = "500 k", "1e+06" = "1 M", "5e+06" = "5 M", "1e+07" = "10 M"), limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("%s", variableName), subtitle = sprintf("%s with effect %g", plotName,q),
         x = "Group size", y = "Sequencing depth",  fill = fillName) # color = "sequencing depth",
  print(heatmap)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_%s_10q%d.pdf", saveName,variableSave,q*10)
    ggsave(filename = path_save, plot = heatmap, height = 5, width = 6, device = cairo_pdf)
    dev.off()
    print(heatmap)}
  rm(heatmap)
}

## FUNCTION for plotting combined mean ROC-plots
# Inputs:   variable = meanROCfinal$m, meanROCfinal$d or meanROCfinal$md 
#           parameterVector = groupSize,sequencingDepth or relations 
#           parameterName = "group size", "sequencing depth" or "relation, trade-off, m*d"
#           parameterSave = "groupsize", "depth", "relation"
#           fillVariable = meanROCfinal$d, meanROCfinal$m or meanROCfinal$plotMD
#           fillName = "Sequencing depth", "Group size" or "Experimental design"
# Outputs:  several plots of combined meanROC-curves which are saved if savePlot==T
plot_combined_meanROCs<-function(plotData, variable, parameterVector, parameterName, parameterSave, fillVariable, fillName, yLim, xLim, strata, strataText, strataName){
  for (i in 1:length(parameterVector)) {
    X=parameterVector[i]
    subtitle=sprintf("Experimental designs with %s %d", parameterName, X, repeats)
    path_save <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%d.pdf", saveName,10*q, parameterSave, X)
    path_save2 <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%d_zoom.pdf", saveName,10*q, parameterSave, X)
    steps= seq(0,1,0.2)
    
    if (all(parameterVector==relations)){
      if (X==3000000){
        Xname = "3 M"
      } else if (X==5000000){
        Xname = "5 M"
      }
      if (strata != 0){
        subtitle=sprintf("Experimental designs with %s %s for %s", Xname, parameterName, strataText)
        path_save <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%d_%s_strata%d.pdf", saveName,10*q, parameterSave, X, strataName, strata)
        path_save2 <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%d_%s_strata%d_zoom.pdf", saveName,10*q, parameterSave, X, strataName, strata)
      } else {
        subtitle=sprintf("Experimental designs with %s %s", Xname, parameterName)
      }
    } else if (all(parameterVector==sequencingDepth)){
      dD=sequencingDepthName[i]
      subtitle=sprintf("Experimental designs with %s %s", parameterName, dD)
      path_save <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%s.pdf", saveName,10*q, parameterSave, dD)
      path_save2 <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%s_zoom.pdf", saveName,10*q, parameterSave, dD)
    }
    
    if (fillName=="Sequencing depth"){
      labelNames = c("10000" = "10 k", "1e+05" = "100 k", "5e+05" = "500 k", "1e+06" = "1 M", "5e+06" = "5 M", "1e+07" = "10 M")
    } else {
      labelNames = waiver()#levels(fillVariable)
    }
      

    if (xLim!=1){
      path_save <- path_save2
      steps= seq(0,xLim,0.005)
    }
    
    combinedPlot<-ggplot(data=plotData[variable==X,], 
                         aes(x=FPR, y=meanTPR, fill=fillVariable[variable==X])) +  #theme_minimal() + 
      geom_ribbon(aes(ymin=(min), ymax=(max),fill = fillVariable[variable==X]), alpha=0.2) +
      geom_line(aes(color = fillVariable[variable==X])) +
      labs(title=sprintf("ROC curves for %s  with effect %g", plotName, q), 
           subtitle = subtitle, x = "False Positive Rate", y = "True Positive Rate",  
           color = fillName, fill = fillName) +
      ylim(0, yLim) + 
      coord_cartesian(xlim=c(0,xLim)) + scale_x_continuous( breaks = steps)+
      scale_fill_viridis_d(begin = 0, end = 0.85, labels = labelNames) + scale_colour_viridis_d(begin = 0, end = 0.85, labels = labelNames)
    print(combinedPlot)
    
    if(savePlot == TRUE){
      ggsave(filename = path_save, plot = combinedPlot, height = 5, width = 6.5)
      dev.off()
      print(combinedPlot)}
    rm(combinedPlot, X)
  }
}
