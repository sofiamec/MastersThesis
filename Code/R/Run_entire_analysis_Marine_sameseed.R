
# Script for running all scripts, performing the entire analysis of datasets with different designs

#===================================================================================================================================
## Loading Libraries:

library(plyr)
library(viridis)
library(DESeq2)
library(ggplot2)
library(pracma)

####################################################################################################################################
#===================================================================================================================================
#                      Only change these parameters for different results! 
#                            (the rest of the code should adjust)
#===================================================================================================================================
## Selecting parameters and data:
onTerra = T                                                # use T if running analysis on Terra (large scale settings applied)
saveName = "Marine"     # "Gut2" or "Marine"                  # this will in turn load the correct data
f = 0.10                                                    # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
numberOfStrata = 3                                          # sets the number of groups for dividing gene abundance and variability
extraDesigns=F                                              # use T if extra designs are added

# Test-settings (CHANGE HERE!)
if (onTerra==F){
  repeats = 3                                               # sets the number of runs for each case (experimental design and q)
  savePlot = F                                              # use T when plots should be saved (for many repeats)
  loadData = F                                              # use T if it is a rerun of existing results
  effectsizes=3#c(1.5,3)                                             # q = Fold-change for downsampling
  groupSize<-c(3,5)#c(3,10)#,10,30,50)                                            # m = Number of samples in each group (total nr samples = 2*m)
  sequencingDepth<-c(10000, 100000)#,1000000)#,500000,1000000,5000000,10000000)      # d = Desired sequencing depth per sample
  sequencingDepthName<-c("10k", "100k")#,"1M")#, "500k","1M","5M","10M")          # dD = Displayed names for sequensing depths
}

# Real settings
if (onTerra==T){
  repeats = 10                                              # sets the number of runs for each case (experimental design and q)
  savePlot = T                                              # use T when plots should be saved (for many repeats)
  loadData = F                                              # use T if it is a rerun of existing results
  effectsizes=c(1.5,3)                                      # q = Fold-change for downsampling
  groupSize<-c(3,5,10,30,50)                                # m = Number of samples in each group (total nr samples = 2*m)
  # sequensing depths are set later depending on dataset    # d and dD = Desired sequencing depths and how it should be displayed
}

# Extra-settings
extraL=0                                                    # unless extraDesigns are added, the length of added designs is 0                               
if (extraDesigns==T){
  # The 3 following must have equal lengths!                # Combined they give more results for trade-off curves. Here with m*d = 3M, 3M and 5M respectively
  extraSeqDepth=c(200000,500000,250000)      
  extraSeqDepthName=c("200k","500k","250k")
  extraGroups=c(15,6,20)
  extraL<-length(extraGroups)
}
#===================================================================================================================================
#===================================================================================================================================
####################################################################################################################################


#===================================================================================================================================
# FUNCTIONS SECTION
#===================================================================================================================================

## FUNCTION to Remove low counts 
# For a given dataset, this function removes genes with low counts (>75 % or an average count <3).
# input: Data = the data to remove genes from
# output: 
remove_low_counts=function(Data){
  a=rowSums(Data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(Data)) {
    b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
    r[i]<-a[i]+b[i]
  }
  FilteredData=Data[r==0,]
  return(FilteredData)
}

## FUNCTION for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

## FUNCTION that assigns factor levels for each gene in a dataset according to variability and abundance
# The created strata will have approximately the same number of genes each  
# Input:    Data = Gut2 or Marine
#           n = the number of strata/levels used when dividing the data
# Output:   DataStrata = a matrix/dataframe? of genes in the given dataset along with their respective stratas (as factors?)
create_strata<-function(Data,numberOfStrata){
  DataStrata<-data.frame(rownames(Data), rowSums(Data),rowVars(as.matrix(Data)))
  n=numberOfStrata
  # create a vector for assigning smalles to highest strata
  repVector<-as.integer(cumsum(c(rep(nrow(Data)/n,n))))
  repVector<-repVector-c(0,repVector)[-(length(repVector)+1)]
  strataVector<-rep(1:n,repVector)
  
  # sort by abundance and assign abundance-strata
  DataStrata<-DataStrata[order(DataStrata[,2]),] 
  DataStrata<-cbind(DataStrata,strataVector)
  
  # sort by variability and assign variability-strata
  DataStrata<-DataStrata[order(DataStrata[,3]),]
  DataStrata<-cbind(DataStrata,strataVector)
  
  # rename columns, factorise strata and sort genes by original order
  colnames(DataStrata)<-c("GeneName", "GeneAbundance", "GeneVariability", "AbundanceStrata", "VariabilityStrata")
  DataStrata$AbundanceStrata<-as.factor(DataStrata$AbundanceStrata)
  DataStrata$VariabilityStrata<-as.factor(DataStrata$VariabilityStrata)
  DataStrata<-DataStrata[rownames(Data),]
  
  rm(n, numberOfStrata, repVector,strataVector)
  return(DataStrata)
}

#===================================================================================================================================
# SET-UP SECTION
#===================================================================================================================================
## Loading original data
#  filtering out samples not meeting expDesign requirements
#  filtering out genes with too low counts

if(saveName == "Gut2"){
  plotName = "Human Gut II"
  Gut2Original <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
  Gut2Intermediate = Gut2Original[,colSums(Gut2Original)>=5000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Gut2 <- remove_low_counts(Gut2Intermediate)
  Data = Gut2
  
  boldvalue2="0"
  relations<-c(3000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000)
    sequencingDepthName<-c("10k","100k","500k","1M", "5M")
  }
  
  rm(Gut2, Gut2Original, Gut2Intermediate)        # remove original and intermediate datasets
  
} else if(saveName == "Marine"){
  plotName = "Marine"
  MarineOriginal <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)
  MarineIntermediate = MarineOriginal[,colSums(MarineOriginal)>=10000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Marine <- remove_low_counts(MarineIntermediate)
  Data = Marine
  
  boldvalue2="5e+07"                                              
  relations<-c(3000000,5000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000,10000000)
    sequencingDepthName<-c("10k","100k","500k","1M", "5M", "10M")
  }
  
  rm(Marine, MarineOriginal, MarineIntermediate)  # remove original and intermediate datasets
}
DataStrata<-create_strata(Data,numberOfStrata)

#===================================================================================================================================
## Setting upp the right environment and naming designs

AllSaveDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 
AllPlotDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 

for (effect in 1:length(effectsizes)) {           # looping over q
  q=effectsizes[effect]
  
  # Creating the standard designs
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      
      AllSaveDesigns[effect,group,seq] <- sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
      AllPlotDesigns[effect,group,seq] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
        
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group,seq]))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group,seq])), recursive = T)
      }
        
      if (!dir.exists(sprintf("../../Result/%s/IntermediatePlots", saveName))){
        dir.create(file.path("../../Result", sprintf("%s", saveName), "/IntermediatePlots"), recursive = T)
      }
    }
  }
  
  # Creating 3 extra designs                              
  if (extraDesigns==T){
  for (i in 1:extraL) {              # looping over extra designs with fixed m and d
    m=extraGroups[i]
    d=extraSeqDepth[i]
    dD=extraSeqDepthName[i]
    
    AllSaveDesigns[effect,group+i,seq+i] <- sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
    AllPlotDesigns[effect,group+i,seq+i] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
    
    # Create folder for certain case if it doesn't exist
    if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group+i,seq+i]))){
      dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group+i,seq+i])), recursive = T)
    }
    
    if (!dir.exists(sprintf("../../Result/%s/IntermediatePlots", saveName))){
      dir.create(file.path("../../Result", sprintf("%s", saveName), "/IntermediatePlots"), recursive = T)
    }
  }}
}


#===================================================================================================================================
# ANALYSIS SECTION
#===================================================================================================================================
## Run the entire analysis for all setups of q, m and d:

set.seed(100)
for (effect in 1:length(effectsizes)) {           # looping over q
  q=effectsizes[effect]
  
  # Creating empty final results-matrices
  meanAUCfinal = data.frame()
  meanROCfinal = data.frame()  
  
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    cat(sprintf("================================== m=%d =======================================\n", m))
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      cat(sprintf("================================== d=%s =====================================\n", dD))
      
      # Creating empty initial result-matrices
      AUC = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),rep = numeric(0))
      ROC = data.frame(TPR=numeric(0), FPR=numeric(0),rep=numeric(0))
      meanROC = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0))
      
      saveExpDesign <- AllSaveDesigns[effect,group,seq]
      plotExpDesign <- AllPlotDesigns[effect,group,seq]
      
      set.seed(100)
      for (run in 1:repeats){
        cat(sprintf("Repeat %d\n", run))
        
        # Run the code for resampling and downsampling, or load already resampled and downsampled dataset
        if (loadData==T){
          DownSampledData<-read.csv(file=sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
          DAGs<-read.csv(file=sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run),header = T,row.names = 1)
          
        } else {      
          source("Resample_datasets_sameseed.R")
        }
        
        # Run the code for analysing DAGs
        source("Analysis_of_DAGs.R")
      
        # Save the results 
        AUC[run,] <- AUCs
        ROC <- rbind(ROC,ROCs)
        meanROC<-rbind(meanROC, meanROCs)
        
        rm(AUCs,ROCs,meanROCs)
      }
      
      ROC$rep<-as.factor(ROC$rep)
      
      # Plot individual ROC-plots
      ROCplot <- ggplot(data=ROC, aes(x=FPR, y=TPR, group=rep)) +  geom_line(aes(color=rep)) + 
        theme(plot.title = element_text(hjust = 0.5)) +  theme_minimal() + 
        scale_color_viridis(begin = 0, end = 0.85, discrete=TRUE) +
        labs(title=sprintf("ROC-curves for %s with effect %g", plotName, q), 
             subtitle = sprintf("Experimental design: %s", plotExpDesign),
             colour="repeats", x = "False Positive Rate", y = "True Positive Rate") +
        xlim(0, 1)+  ylim(0, 1)
      
      print(ROCplot)
      
      if(savePlot == TRUE){
        path_save <-  sprintf("../../Result/%s/IntermediatePlots/individualROCs_rep%s.pdf", saveName, saveExpDesign, rep)
        ggsave(filename = path_save, plot = ROCplot, height = 5, width = 6)
        dev.off()
        print(ROCplot)
      }
      
      colnames(meanROC)[3]<-"meanTPR"
      meanROC2<-ddply(meanROC, "FPR", summarise,
                      N    = length(meanTPR),
                      mean = mean(meanTPR),
                      min  = min(meanTPR),
                      max  = max(meanTPR))
      colnames(meanROC2)[3]<-"meanTPR"
      
      # Plot mean RoC-curves for certain experimental design
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
      }
      meanAUCfinal<-rbind(meanAUCfinal,data.frame(t(colMeans(AUC)[1:5]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      
      rm(ROCplot, meanROC, meanROCplot, dD)
      rm(ROC, AUC,  meanROC2, run)
    }
  }
  
  # Repeating analysis for extra designs #================================================#
  if(extraDesigns == TRUE){
  for (i in 1:extraL){ # looping fixed m and d
    m=extraGroups[i]
    cat(sprintf("================================== m=%d =======================================\n", m))
    d=extraSeqDepth[i]
    dD=extraSeqDepthName[i]
    cat(sprintf("================================== d=%s =====================================\n", dD))
      
    # Creating empty initial result-matrices
    AUC = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),rep = numeric(0))
    ROC = data.frame(TPR=numeric(0), FPR=numeric(0),rep=numeric(0))
    meanROC = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0))
      
    saveExpDesign <- AllSaveDesigns[effect,group+i,seq+i]
    plotExpDesign <- AllPlotDesigns[effect,group+i,seq+i]
      
    for (run in 1:repeats){
      cat(sprintf("Repeat %d\n", run))
        
      # Run the code for resampling and downsampling, 
      # or load already resampled and downsampled dataset
      if (loadData==T){
        DownSampledData<-read.csv(file=sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
        DAGs<-read.csv(file=sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run),header = T,row.names = 1)
        
      } else {      
        source("Resample_datasets.R")
      }
      
      # Run the code for analysing DAGs
      source("Analysis_of_DAGs.R")
        
      # Save the results 
      AUC[run,] <- AUCs
      ROC <- rbind(ROC,ROCs)
      meanROC<-rbind(meanROC, meanROCs)
        
      rm(AUCs,ROCs,meanROCs)
    }
    
    ROC$rep<-as.factor(ROC$rep)
    
    # Plot individual ROC-plots
    ROCplot <- ggplot(data=ROC, aes(x=FPR, y=TPR, group=rep)) +  geom_line(aes(color=rep)) + 
      theme(plot.title = element_text(hjust = 0.5)) +  theme_minimal() + 
      scale_color_viridis(begin = 0, end = 0.85, discrete=TRUE) +
      labs(title=sprintf("ROC-curves for %s with effect %g", plotName, q), 
           subtitle = sprintf("Experimental design: %s", plotExpDesign),
           colour="repeats", x = "False Positive Rate", y = "True Positive Rate") +
      xlim(0, 1)+  ylim(0, 1)
    
    print(ROCplot)
      
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/IntermediatePlots/individualROCs_rep%s.pdf", saveName, saveExpDesign, rep)
      ggsave(filename = path_save, plot = ROCplot, height = 5, width = 6)
      dev.off()
      print(ROCplot)
    }
      
    colnames(meanROC)[3]<-"meanTPR"
    meanROC2<-ddply(meanROC, "FPR", summarise,
                    N    = length(meanTPR),
                    mean = mean(meanTPR),
                    min  = min(meanTPR),
                    max  = max(meanTPR))
    colnames(meanROC2)[3]<-"meanTPR"
      
    # Plot mean RoC-curves for certain experimental design
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
    }
    meanAUCfinal<-rbind(meanAUCfinal,data.frame(t(colMeans(AUC)[1:5]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
    meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      
    rm(ROCplot, meanROC, meanROCplot, dD)
    rm(ROC, AUC,  meanROC2, run)
  }}
  
  #===================================================================================================================================
  ### Summarising results:
  
  colnames(meanAUCfinal)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "d", "m" ,"md","plotMD")
  colnames(meanROCfinal)<-c("FPR", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
  
  meanAUCfinal$md[meanAUCfinal$md=="5e+06"]<-"bold"
  meanAUCfinal$md[meanAUCfinal$md==boldvalue2]<-"bold"
  meanAUCfinal$md[meanAUCfinal$md!="bold"]<-"plain"
  if (extraDesigns==T){
    HeatmapData<-head(meanAUCfinal,-extraL) 
  } else if (extraDesigns==F){
    HeatmapData<-meanAUCfinal 
  }
  HeatmapData$d<-as.factor(HeatmapData$d)
  HeatmapData$m<-as.factor(HeatmapData$m)
  HeatmapData$md<-as.factor(HeatmapData$md)
  
  meanROCfinal$d<-as.factor(meanROCfinal$d)
  meanROCfinal$m<-as.factor(meanROCfinal$m)
  meanROCfinal$plotMD = factor(meanROCfinal$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)
  
  # Save tables:
  write.csv(meanAUCfinal, file=sprintf("../../Result/%s/AUC_10q%d.csv", saveName,10*q))
  
  # heatmaps for AUC and TPR at FPR 0.01, 0.05 and 1
  heatmapAUCtot <- ggplot(HeatmapData, aes(x=HeatmapData$m, y=HeatmapData$d, fill=AUCtot)) +
    geom_tile(aes(fill = AUCtot)) + geom_text(aes(label = round2(AUCtot, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean total AUC-values for %s with effect %g", plotName,q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
  print(heatmapAUCtot)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_AUCtot_10q%d.pdf", saveName,q*10)
    ggsave(filename = path_save, plot = heatmapAUCtot, height = 5, width = 6)
    dev.off()
    print(heatmapAUCtot)}
  
  heatmapAUC5 <- ggplot(HeatmapData, aes(x=HeatmapData$m, y=HeatmapData$d, fill=AUC5)) +
    geom_tile(aes(fill = AUC5)) + geom_text(aes(label = round2(AUC5, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean AUC-values at FPR 0.05 for %s with effect %g", plotName,q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
  print(heatmapAUC5)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_AUC5_10q%d.pdf", saveName,10*q)
    ggsave(filename = path_save, plot = heatmapAUC5, height = 5, width = 6)
    dev.off()
    print(heatmapAUC5)}
  
  heatmapAUC1 <- ggplot(HeatmapData, aes(x=HeatmapData$m, y=HeatmapData$d, fill=AUC1)) +
    geom_tile(aes(fill = AUC1)) + geom_text(aes(label = round2(AUC1, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean AUC-values at FPR 0.01 for %s  with effect %g", plotName, q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
  print(heatmapAUC1)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_AUC1_10q%d.pdf", saveName,q*10)
    ggsave(filename = path_save, plot = heatmapAUC1, height = 5, width = 6)
    dev.off()
    print(heatmapAUC1)}
  
  heatmapTPR5 <- ggplot(HeatmapData, aes(x=HeatmapData$m, y=HeatmapData$d, fill=TPR5)) +
    geom_tile(aes(fill = TPR5)) + geom_text(aes(label = round2(TPR5, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean TPR-values at FPR 0.05 for %s with effect %g", plotName, q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "TPR-values") 
  print(heatmapTPR5)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_TPR5_10q%d.pdf", saveName,q*10)
    ggsave(filename = path_save, plot = heatmapTPR5, height = 5, width = 6)
    dev.off()
    print(heatmapTPR5)}
  
  heatmapTPR1 <- ggplot(HeatmapData, aes(x=HeatmapData$m, y=HeatmapData$d, fill=TPR1)) +
    geom_tile(aes(fill = TPR1)) + geom_text(aes(label = round2(TPR1, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean TPR-values at FPR 0.01 for %s with effect %g", plotName, q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "TPR-values") 
  print(heatmapTPR1)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_TPR1_10q%d.pdf", saveName,q*10)
    ggsave(filename = path_save, plot = heatmapTPR1, height = 5, width = 6)
    dev.off()
    print(heatmapTPR1)}
  
  
  ### Plot mean RoC-curves for all experimental designs
  
  # mean plots with set groupsize
  for (group in 1:length(groupSize)){
    M=groupSize[group]
    
    meanROCplotgroup <- ggplot(data=meanROCfinal[meanROCfinal$m==M,], aes(x=FPR, y=meanTPR, fill=d)) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(min), ymax=(max),fill = d), alpha=0.2) +
      geom_line(aes(color = d)) +
      labs(title=sprintf("Mean ROC-curves for %s  with effect %g", plotName, q), 
           subtitle = sprintf("Experimental designs with group size %d     (%s repeats each)", M, repeats),
           x = "False Positive Rate", y = "True Positive Rate",  color = "sequensing depth", fill = "sequensing depth") +
      ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_fill_viridis_d(begin = 0, end = 0.85) +
      scale_colour_viridis_d(begin = 0, end = 0.85)
    
    print(meanROCplotgroup)
    
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/meanROC_10q%d_groupsize_%d.pdf", saveName,10*q, M)
      ggsave(filename = path_save, plot = meanROCplotgroup, height = 5, width = 6)
      dev.off()
      print(meanROCplotgroup)}
    
    rm(meanROCplotgroup, M)
  }
  
  # mean plots with set sequencing depth
  for (seq in 1:length(sequencingDepth)) {
    D=sequencingDepth[seq]
    dD=sequencingDepthName[seq]
    
    meanROCplotdepth <- ggplot(data=meanROCfinal[meanROCfinal$d==D,], aes(x=FPR, y=meanTPR, fill=m)) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(min), ymax=(max),fill = m), alpha=0.2) +
      geom_line(aes(color = m)) +
      labs(title=sprintf("Mean ROC-curves for %s with effect %g", plotName, q), 
           subtitle = sprintf("Experimental designs with sequencing depth %s     (%s repeats each)", dD, repeats),
           x = "False Positive Rate", y = "True Positive Rate",  color = "group size", fill = "group size") +
      ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_fill_viridis_d(begin = 0, end = 0.85) +
      scale_colour_viridis_d(begin = 0, end = 0.85)
    
    print(meanROCplotdepth)
    
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/meanROC_10q%d_depth_%d.pdf", saveName,10*q,D)
      ggsave(filename = path_save, plot = meanROCplotdepth, height = 5, width = 6)
      dev.off()
      print(meanROCplotdepth)}
    
    rm(meanROCplotdepth, D,dD)
  }
  
  # mean plots with set groupsize and depth relation (m*d)
  for (relation in 1:length(relations)){
    MD=relations[relation]
    
    meanROCplotrelation <- ggplot(data=meanROCfinal[meanROCfinal$md==MD,], aes(x=FPR, y=meanTPR, fill=plotMD)) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(min), ymax=(max),fill = plotMD), alpha=0.2) +
      geom_line(aes(color = plotMD)) +
      labs(title=sprintf("Mean ROC-curves for %s  with effect %g", plotName, q), 
           subtitle = sprintf("Experimental designs with relation/trade-off/m*d=%d     (%s repeats each)", MD, repeats),
           x = "False Positive Rate", y = "True Positive Rate",  color = "Experimental design", fill = "Experimental design") +
      ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_fill_viridis_d(begin = 0, end = 0.85) +
      scale_colour_viridis_d(begin = 0, end = 0.85)
    
    print(meanROCplotrelation)
    
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/meanROC_10q%d_relation_%d.pdf", saveName,10*q, MD)
      ggsave(filename = path_save, plot = meanROCplotrelation, height = 5, width = 6)
      dev.off()
      print(meanROCplotrelation)}
    
    rm(meanROCplotrelation, MD)
  }
  
  rm(group,seq, relation)

}

rm(repeats)
