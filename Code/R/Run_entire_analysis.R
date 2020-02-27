
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
onTerra = F                                                 # use T if running analysis on Terra (large scale settings applied)
saveName = "Gut2"     # "Gut2" or "Marine"                  # this will in turn load the correct data
f = 0.10                                                    # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
numberOfStrata = 3                                          # sets the number of groups for dividing gene abundance and variability
extraDesigns=F                                              # use T if extra designs are added

# Test-settings (CHANGE HERE!)
if (onTerra==F){
  repeats = 2                                               # sets the number of runs for each case (experimental design and q)
  savePlot = F                                              # use T when plots should be saved (for many repeats)
  loadData = T                                              # use T if it is a rerun of existing results
  effectsizes=3                                             # q = Fold-change for downsampling
  groupSize<-c(3,5)#,10,30,50)                                            # m = Number of samples in each group (total nr samples = 2*m)
  sequencingDepth<-c(10000,100000)#,500000,1000000,5000000,10000000)      # d = Desired sequencing depth per sample
  sequencingDepthName<-c("10k","100k")#, "500k","1M","5M","10M")          # dD = Displayed names for sequensing depths
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

## FUNCTION for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

## FUNCTION for plotting AUC- and TPR-heatmaps
# Inputs:   variable = HeatmapData$AUCtot,HeatmapData$TPR1 etc
#           variableName = "total AUC-values", "TPR-values at FPR 0.01" etc
#           fillName = "AUC-values" or "TPR-values"
#           variableSave = "AUCtot", "TPR1" etc
# Outputs:  a plotted heatmap which is saved if savePlot==T
plot_heatmaps<-function(variable,variableName, fillName, variableSave){
  heatmap <- ggplot(HeatmapData, aes(x=m, y=d, fill=variable)) +
    geom_tile(aes(fill = variable)) + geom_text(aes(label = round2(variable, 2), fontface=md)) +
    scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
    scale_y_discrete(limits = rev(levels(as.factor(HeatmapData$d)))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
    labs(title=sprintf("Mean %s for %s with effect %g", variableName, plotName,q), 
         x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = fillName) 
  print(heatmap)
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/heatmap_%s_10q%d.pdf", variableSave, saveName,q*10)
    ggsave(filename = path_save, plot = heatmap, height = 5, width = 6)
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
plot_combined_meanROCs<-function(variable, parameterVector, parameterName, parameterSave, fillVariable, fillName){
  for (i in 1:length(parameterVector)) {
    X=parameterVector[i]
    subtitle=sprintf("Experimental designs with %s %d     (%s repeats each)", parameterName, X, repeats)
    path_save <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%d.pdf", saveName,10*q, parameterSave, X)
    if (all(parameterVector==sequencingDepth)){
      dD=sequencingDepthName[i]
      subtitle=sprintf("Experimental designs with %s %s     (%s repeats each)", parameterName, dD, repeats)
      path_save <-  sprintf("../../Result/%s/meanROC_10q%d_%s_%s.pdf", saveName,10*q, parameterSave, dD)
    }
    combinedPlot<-ggplot(data=meanROCfinal[variable==X,], aes(x=FPR, y=meanTPR, fill=fillVariable[variable==X])) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(min), ymax=(max),fill = fillVariable[variable==X]), alpha=0.2) +
      geom_line(aes(color = fillVariable[variable==X])) +
      labs(title=sprintf("Mean ROC-curves for %s  with effect %g", plotName, q), 
           subtitle = subtitle, x = "False Positive Rate", y = "True Positive Rate",  
           color = fillName, fill = fillName) +
      ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      scale_fill_viridis_d(begin = 0, end = 0.85) + scale_colour_viridis_d(begin = 0, end = 0.85)
    print(combinedPlot)
    if(savePlot == TRUE){
      ggsave(filename = path_save, plot = combinedPlot, height = 5, width = 6)
      dev.off()
      print(combinedPlot)}
    rm(combinedPlot, X)
  }
}


#===================================================================================================================================
# SET-UP SECTION
#===================================================================================================================================

source("Setup_designs.R")

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
      
      for (run in 1:repeats){
        cat(sprintf("Repeat %d\n", run))
        
        # Run the code for resampling and downsampling, or load already resampled and downsampled dataset
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
        
        # STRATA
        #ROCsAbundance
        #AUCsAbundance 
        #meanROCsAbundance
        
        #ROCsVariability
        #AUCsVariability
        #meanROCsVariability 
        
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
  
  ### Plotting heatmaps for AUC- and TPR-values
  plot_heatmaps(HeatmapData$AUC1, "AUC-values at FPR 0.01", "AUC-values", "AUC1")
  plot_heatmaps(HeatmapData$AUC5, "AUC-values at FPR 0.05", "AUC-values", "AUC5")
  plot_heatmaps(HeatmapData$AUCtot, "total AUC-values", "AUC-values", "AUCtot")
  plot_heatmaps(HeatmapData$TPR1, "TPR-values at FPR 0.01", "TPR-values", "TPR1")
  plot_heatmaps(HeatmapData$TPR5, "TPR-values at FPR 0.05", "TPR-values", "TPR5")
  
  ### Plot mean RoC-curves for all experimental designs
  
  # mean plots with set groupsize
  plot_combined_meanROCs(meanROCfinal$m, groupSize, "group size", "groupsize", meanROCfinal$d, "Sequencing depth")
  # mean plots with set sequencing depth
  plot_combined_meanROCs(meanROCfinal$d, sequencingDepth, "sequencing depth", "depth", meanROCfinal$m, "Group size")
  # mean plots with set groupsize and depth relation (m*d)
  plot_combined_meanROCs(meanROCfinal$md, relations, "relation/trade-off/m*d", "relation", meanROCfinal$plotMD, "Experimental design")
  
  rm(group,seq, relation)

}

rm(repeats)
