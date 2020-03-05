
# Script for running all scripts, performing the entire analysis of datasets with different designs

## Nu kallar denna på analysis_of_DAGs_strata och REsample_data_strata!

####################################################################################################################################
#===================================================================================================================================
#                      Only change these parameters for different results! 
#                            (the rest of the code should adjust)
#===================================================================================================================================
## Selecting parameters and data:
onTerra = T                                                 # use T if running analysis on Terra (large scale settings applied)
saveName = "Gut2"     # "Gut2" or "Marine"                  # this will in turn load the correct data
f = 0.10                                                    # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
runStrata = T
extraDesigns=T                                              # use T if extra designs are added

# Test-settings (CHANGE HERE!)
if (onTerra==F){
  repeats = 2                                               # sets the number of runs for each case (experimental design and q)
  savePlot = T                                              # use T when plots should be saved (for many repeats)
  loadData = F                                              # use T if it is a rerun of existing results
  effectsizes=1.5                                             # q = Fold-change for downsampling
  groupSize<-c(3)#,30)#,10,30,50)                                            # m = Number of samples in each group (total nr samples = 2*m)
  sequencingDepth<-c(10000,100000,500000)#,500000,1000000,5000000)      # d = Desired sequencing depth per sample
  sequencingDepthName<-c("10k","100k","500k")#, "500k","1M","5M")          # dD = Displayed names for sequensing depths
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
 # Strata-settings
if (runStrata==T){
  numberOfStrata = 3                                          # sets the number of groups for dividing gene abundance and variability
  strataClass<-factor(c("low","medium","high"), levels=c("low","medium","high")) #strataClass<-c("low","medium", "high")                      # should correspond to the number of stratas
}
#===================================================================================================================================
#===================================================================================================================================
####################################################################################################################################

#===================================================================================================================================
# FUNCTIONS and LIBRARIES SECTION
#===================================================================================================================================

source("Functions_and_libraries.R")

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
  meanGenesFDR = data.frame()
  if (runStrata==T){
    #meanAUCfinalAb = data.frame()
    #meanAUCfinalV = data.frame()
    meanROCfinalAb = data.frame() 
    meanROCfinalV = data.frame() 
  }
  
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    cat(sprintf("================================== m=%d =======================================\n", m))
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      cat(sprintf("================================== d=%s =====================================\n", dD))
      
      # Creating empty initial result-matrices
      AUC = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0))
      ROC = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0))
      meanROC = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0))
      genesFDR = data.frame()
      if (runStrata==T){
        #AUCAbundance = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
        #AUCVariability = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
        ROCAbundance = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0), strata=numeric(0))
        ROCVariability = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0), strata=numeric(0))
        meanROCAbundance = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), strata=numeric(0))
        meanROCVariability = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), strata=numeric(0))
      }
      
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
        meanROC <- rbind(meanROC, meanROCs)
        genesFDR <- rbind(genesFDR,genesFDRs)
        
        if (runStrata==T){
          #AUCAbundance<-rbind(AUCAbundance, AUCsAbundance)
          #AUCVariability<-rbind(AUCVariability, AUCsVariability)
          ROCAbundance <- rbind(ROCAbundance,ROCsAbundance)
          ROCVariability <- rbind(ROCVariability,ROCsVariability)
          meanROCAbundance <- rbind(meanROCAbundance,meanROCsAbundance)
          meanROCVariability <- rbind(meanROCVariability,meanROCsVariability)
          rm(ROCsAbundance, ROCsVariability, meanROCsAbundance, meanROCsVariability) # rm(AUCsAbundance, AUCsVariability)
        }
        
      rm(AUCs, ROCs, meanROCs)
      }
      
      # plot individual ROCs
      ROC$run<-as.factor(ROC$run)
      individual_ROC_plot(ROC)

      colnames(meanROC)[3]<-"meanTPR"
      meanROC2<-ddply(meanROC, "FPR", summarise,
                      N    = length(meanTPR),
                      mean = mean(meanTPR),
                      min  = min(meanTPR),
                      max  = max(meanTPR))
      colnames(meanROC2)[3]<-"meanTPR"
      
      # Save results for this experimental design
      meanAUCfinal<-rbind(meanAUCfinal,data.frame(t(colMeans(AUC)[1:5]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      meanGenesFDR<-rbind(meanGenesFDR,data.frame(t(colMeans(genesFDR)),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      
      if(runStrata==T){
        ROCAbundance$run<-as.factor(ROCAbundance$run)
        ROCVariability$run<-as.factor(ROCVariability$run)
        ROCAbundance$strata<-as.factor(ROCAbundance$strata)
        ROCVariability$strata<-as.factor(ROCVariability$strata)
        
        meanROCAb2<-ddply(meanROCAbundance, c("FPR", "strata"), summarise,
                          N    = length(meanTPR),
                          mean = mean(meanTPR),
                          min  = min(meanTPR),
                          max  = max(meanTPR))
        
        meanROCV2<-ddply(meanROCVariability, c("FPR", "strata"), summarise,
                         N    = length(meanTPR),
                         mean = mean(meanTPR),
                         min  = min(meanTPR),
                         max  = max(meanTPR))
        
        colnames(meanROCAb2)[4]<-"meanTPR"
        colnames(meanROCV2)[4]<-"meanTPR"
        
        meanROCAb2$strata<-as.factor(meanROCAb2$strata)
        meanROCV2$strata<-as.factor(meanROCV2$strata)
        
        # Plot mean RoC-curves for certain experimental design with Abundance-strata
        meanROCplotAb <- ggplot(data=meanROCAb2, aes(x=FPR, y=meanTPR, group=strata)) + 
          geom_ribbon(aes(ymin=(min), ymax=(max), fill=strataClass[strata]), alpha = 0.2) + 
          geom_line(aes(color=strataClass[strata])) +   theme_minimal() + 
          labs(title=sprintf("Mean ROC-curve for %s with effect %g", plotName,q), 
               subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
               x = "False Positive Rate", y = "True Positive Rate", color="Abundance strata", fill="Abundance strata")+
          ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
          scale_fill_viridis_d(begin = 0.2, end = 0.6) +
          scale_colour_viridis_d(begin = 0.2, end = 0.6)
        print(meanROCplotAb)
        if(savePlot == TRUE){
          path_save <-  sprintf("../../Result/%s/IntermediatePlots/meanROC_%s_AbStrata.pdf", saveName, saveExpDesign)
          ggsave(filename = path_save, plot = meanROCplotAb, height = 5, width = 6)
          dev.off()
          print(meanROCplotAb)
          rm(path_save)
        }
        rm(meanROCplotAb)
        
        # Plot mean RoC-curves for certain experimental design with Variability-strata
        meanROCplotV <- ggplot(data=meanROCV2, aes(x=FPR, y=meanTPR, group=strata)) + 
          geom_ribbon(aes(ymin=(min), ymax=(max), fill=strataClass[strata]), alpha = 0.2) + 
          geom_line(aes(color=strataClass[strata])) + theme_minimal() + 
          labs(title=sprintf("Mean ROC-curve for %s with effect %g", plotName,q), 
               subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
               x = "False Positive Rate", y = "True Positive Rate", color="Variability strata", fill="Variability strata")+
          ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
          scale_fill_viridis_d(begin = 0.2, end = 0.6) +
          scale_colour_viridis_d(begin = 0.2, end = 0.6)
        print(meanROCplotV)
        if(savePlot == TRUE){
          path_save <-  sprintf("../../Result/%s/IntermediatePlots/meanROC_%s_VarStrata.pdf", saveName, saveExpDesign)
          ggsave(filename = path_save, plot = meanROCplotV, height = 5, width = 6)
          dev.off()
          print(meanROCplotV)
          rm(path_save)
        }
        rm(meanROCplotV)
        
        # Save strata-results for this experimental design
        #for (k in 1:numberOfStrata) {
        #  meanAUCfinalAb <-rbind(meanAUCfinalAb,data.frame(t(colMeans(AUCAbundance[AUCAbundance$strata==k,])[-6]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
        #  meanAUCfinalV <- rbind(meanAUCfinalV,data.frame(t(colMeans(AUCVariability[AUCVariability$strata==k,])[-6]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
        #}
        meanROCfinalAb<-rbind(meanROCfinalAb,data.frame(meanROCAb2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
        meanROCfinalV<-rbind(meanROCfinalV,data.frame(meanROCV2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
        rm(meanROCAbundance, meanROCVariability, meanROCAb2, meanROCV2)
      }
      
      rm(meanROC, dD)
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
    genesFDR = data.frame()
    if (runStrata==T){
      ROCAbundance = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0), strata=numeric(0))
      ROCVariability = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0), strata=numeric(0))
      meanROCAbundance = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), strata=numeric(0))
      meanROCVariability = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), strata=numeric(0))
    }
    
    saveExpDesign <- AllSaveDesigns[effect,group+i,seq+i]
    plotExpDesign <- AllPlotDesigns[effect,group+i,seq+i]
      
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
      genesFDR <- rbind(genesFDR,genesFDRs)
      if (runStrata==T){
        ROCAbundance <- rbind(ROCAbundance,ROCsAbundance)
        ROCVariability <- rbind(ROCVariability,ROCsVariability)
        meanROCAbundance <- rbind(meanROCAbundance,meanROCsAbundance)
        meanROCVariability <- rbind(meanROCVariability,meanROCsVariability)
        rm(ROCsAbundance, ROCsVariability, meanROCsAbundance, meanROCsVariability)
      }
      
      rm(AUCs, ROCs, meanROCs, genesFDRs)
    }
    
    # plot individual ROCs
    ROC$run<-as.factor(ROC$run)
    individual_ROC_plot(ROC)
    
    colnames(meanROC)[3]<-"meanTPR"
    meanROC2<-ddply(meanROC, "FPR", summarise,
                    N    = length(meanTPR),
                    mean = mean(meanTPR),
                    min  = min(meanTPR),
                    max  = max(meanTPR))
    colnames(meanROC2)[3]<-"meanTPR"
    
    # Save results for this experimental design
    meanAUCfinal<-rbind(meanAUCfinal,data.frame(t(colMeans(AUC)[1:5]),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
    meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
    meanGenesFDR<-rbind(meanGenesFDR,data.frame(t(colMeans(genesFDR)),d,m,m*d,sprintf("m=%d d=%s",m,dD)))
    
    if (runStrata==T){
      ROCAbundance$run<-as.factor(ROCAbundance$run)
      ROCVariability$run<-as.factor(ROCVariability$run)
      ROCAbundance$strata<-as.factor(ROCAbundance$strata)
      ROCVariability$strata<-as.factor(ROCVariability$strata)
    
      # No individual plots for strata in extra designs
    
      meanROCAb2<-ddply(meanROCAbundance, c("FPR", "strata"), summarise,
                        N    = length(meanTPR),
                        mean = mean(meanTPR),
                        min  = min(meanTPR),
                        max  = max(meanTPR))
      
      meanROCV2<-ddply(meanROCVariability, c("FPR", "strata"), summarise,
                       N    = length(meanTPR),
                       mean = mean(meanTPR),
                       min  = min(meanTPR),
                       max  = max(meanTPR))
      
      colnames(meanROCAb2)[4]<-"meanTPR"
      colnames(meanROCV2)[4]<-"meanTPR"
      
      meanROCAb2$strata<-as.factor(meanROCAb2$strata)
      meanROCV2$strata<-as.factor(meanROCV2$strata)
      
      # Save strata-results for this experimental design
      meanROCfinalAb <- rbind(meanROCfinalAb,data.frame(meanROCAb2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      meanROCfinalV <- rbind(meanROCfinalV,data.frame(meanROCV2,d,m,m*d,sprintf("m=%d d=%s",m,dD)))
      
      rm(ROCAbundance, ROCVariability, meanROCAbundance, meanROCVariability, meanROCAb2, meanROCV2)
    }

    rm(ROC, AUC,  meanROC2, run, meanROC, genesFDR)
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
  
  ## for strata
  if (runStrata==T){
    #colnames(meanAUCfinalAb)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "strata", "d", "m" ,"md","plotMD")
    colnames(meanROCfinalAb)<-c("FPR", "strata", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
    #colnames(meanAUCfinalV)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "strata", "d", "m" ,"md","plotMD")
    colnames(meanROCfinalV)<-c("FPR", "strata", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
    
    meanROCfinalAb$d<-as.factor(meanROCfinalAb$d)
    meanROCfinalAb$m<-as.factor(meanROCfinalAb$m)
    meanROCfinalAb$plotMD = factor(meanROCfinalAb$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)
    
    meanROCfinalV$d<-as.factor(meanROCfinalV$d)
    meanROCfinalV$m<-as.factor(meanROCfinalV$m)
    meanROCfinalV$plotMD = factor(meanROCfinalV$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)
  }

  # Save tables:
  write.csv(meanAUCfinal, file=sprintf("../../Result/%s/AUC_10q%d.csv", saveName,10*q))
  write.csv(meanGenesFDR, file=sprintf("../../Result/%s/GenesFDR_10q%d.csv", saveName,10*q))
  
  ### Plotting heatmaps for AUC- and TPR-values
  plot_heatmaps(HeatmapData$AUC1, "AUC-values at FPR 0.01", "AUC-values", "AUC1")
  plot_heatmaps(HeatmapData$AUC5, "AUC-values at FPR 0.05", "AUC-values", "AUC5")
  plot_heatmaps(HeatmapData$AUCtot, "total AUC-values", "AUC-values", "AUCtot")
  plot_heatmaps(HeatmapData$TPR1, "TPR-values at FPR 0.01", "TPR-values", "TPR1")
  plot_heatmaps(HeatmapData$TPR5, "TPR-values at FPR 0.05", "TPR-values", "TPR5")
  
  ### Plot mean RoC-curves for all experimental designs
  # mean plots with set groupsize
  plot_combined_meanROCs(meanROCfinal, meanROCfinal$m, groupSize, "group size", "groupsize", meanROCfinal$d, "Sequencing depth", 1, 1, 0)
  # mean plots with set sequencing depth
  plot_combined_meanROCs(meanROCfinal, meanROCfinal$d, sequencingDepth, "sequencing depth", "depth", meanROCfinal$m, "Group size", 1, 1, 0)
  # mean plots with set groupsize and depth relation (m*d)
  plot_combined_meanROCs(meanROCfinal, meanROCfinal$md, relations, "relation/trade-off/m*d", "relation", meanROCfinal$plotMD, "Experimental design", 1, 1, 0)
  plot_combined_meanROCs(meanROCfinal, meanROCfinal$md, relations, "relation/trade-off/m*d", "relation", meanROCfinal$plotMD, "Experimental design", 1, 0.01, 0)
  
  
  ### Plot mean ROC-curves for strata
  if (runStrata==T){
  # mean plots with set groupsize and depth relation (m*d)
    for (strata in 1:numberOfStrata){
      class<-strataClass[strata]
      plotData=meanROCfinalAb[meanROCfinalAb$strata==strata,]
      plot_combined_meanROCs(plotData, plotData$md, relations, "relation/trade-off/m*d", "relation", plotData$plotMD, "Experimental design", 1, 1, strata, sprintf("genes with %s abundance", class))
      plotData=meanROCfinalV[meanROCfinalV$strata==strata,]
      plot_combined_meanROCs(plotData, plotData$md, relations, "relation/trade-off/m*d", "relation", plotData$plotMD, "Experimental design", 1, 1, strata, sprintf("genes with %s variability", class))
    }
    rm(strata, plotData, class)
  }
  
  rm(group,seq, m, d, dD)
}

rm(repeats, effect, q, extraL, f, relations, boldvalue2, AllPlotDesigns, AllSaveDesigns, plotExpDesign, saveExpDesign)
rm(compute_low_counts, Compute_ROC_AUC, DESeq2_analysis, individual_ROC_plot,introducing_DAGs, plot_combined_meanROCs,plot_heatmaps, resample, round2)

