
# Script for running all scripts, performing the entire analysis of datasets with different designs

#===================================================================================================================================
# ANALYSIS SECTION
#===================================================================================================================================
## Run the entire analysis for all setups of q, m and d:

set.seed(100)

for (effect in 1:length(effectsizes)) {           # looping over q
  q=effectsizes[effect]
  
  # Creating empty final results-matrices
  medianAUCfinal = data.frame()
  meanROCfinal = data.frame() 
  medianGenesFDR = data.frame()
  if (runStrata==T){
    medianAUCfinalAb = data.frame()
    medianAUCfinalV = data.frame()
    meanROCfinalAb = data.frame() 
    meanROCfinalV = data.frame() 
  }
  
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    cat(sprintf("================================== m=%d =======================================\n", m))
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      dS=sequencingSaveName[seq]
      cat(sprintf("================================== d=%s =====================================\n", dD))
      
      # Creating empty initial result-matrices
      AUC = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0))
      ROC = data.frame(TPR=numeric(0), FPR=numeric(0),run=numeric(0))
      meanROC = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0))
      genesFDR = data.frame()
      if (runStrata==T){
        AUCAbundance = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
        AUCVariability = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
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
          source("Resample_datasets.R", local = TRUE)
        }
        
        # Run the code for analysing DAGs
        source("Analysis_of_DAGs.R", local = TRUE)
      
        # Save the results 
        AUC[run,] <- AUCs
        ROC <- rbind(ROC,ROCs)
        meanROC <- rbind(meanROC, meanROCs)
        genesFDR <- rbind(genesFDR,genesFDRs)
        
        if (runStrata==T){
          AUCAbundance<-rbind(AUCAbundance, AUCsAbundance)
          AUCVariability<-rbind(AUCVariability, AUCsVariability)
          ROCAbundance <- rbind(ROCAbundance,ROCsAbundance)
          ROCVariability <- rbind(ROCVariability,ROCsVariability)
          meanROCAbundance <- rbind(meanROCAbundance,meanROCsAbundance)
          meanROCVariability <- rbind(meanROCVariability,meanROCsVariability)
          rm(ROCsAbundance, ROCsVariability, meanROCsAbundance, meanROCsVariability, AUCsAbundance, AUCsVariability)
        }
        
      rm(AUCs, ROCs, meanROCs)
      }
      
      # plot individual ROCs
      #ROC$run<-as.factor(ROC$run)
      #individual_ROC_plot(ROC)

      colnames(meanROC)[3]<-"meanTPR"
      meanROC2<-ddply(meanROC, "FPR", summarise,
                      N    = length(meanTPR),
                      mean = mean(meanTPR),
                      min  = min(meanTPR),
                      max  = max(meanTPR))
      colnames(meanROC2)[3]<-"meanTPR"
      
      # Save results for this experimental design
      medianAUCfinal<-rbind(medianAUCfinal,data.frame(t(colMedians(as.matrix(AUC)[,1:5],na.rm = T)),t(colSums(!is.na(AUC[1:5]))),d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
      meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
      medianGenesFDR<-rbind(medianGenesFDR,data.frame(t(colMedians(as.matrix(genesFDR),na.rm = T)), t(colSums(!is.na(genesFDR))),sprintf("m=%d, d=%s",m,dD)))
      
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
        #meanROCplotAb <- ggplot(data=meanROCAb2, aes(x=FPR, y=meanTPR, group=strata)) + 
        #  geom_ribbon(aes(ymin=(min), ymax=(max), fill=strataClass[strata]), alpha = 0.2) + 
        #  geom_line(aes(color=strataClass[strata])) + #theme_minimal() + 
        #  labs(title=sprintf("Mean ROC-curve for %s with effect %g", plotName,q), 
        #       subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
        #       x = "False Positive Rate", y = "True Positive Rate", color="Abundance strata", fill="Abundance strata")+
        #  ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
        #  scale_fill_viridis_d(begin = 0.2, end = 0.6) +
        #  scale_colour_viridis_d(begin = 0.2, end = 0.6)
        #print(meanROCplotAb)
        #if(savePlot == TRUE){
        #  path_save <-  sprintf("../../Result/%s/IntermediatePlots/meanROC_%s_AbStrata.pdf", saveName, saveExpDesign)
        #  ggsave(filename = path_save, plot = meanROCplotAb, height = 5, width = 6)
        #  dev.off()
        #  print(meanROCplotAb)
        #  rm(path_save)
        #}
        #rm(meanROCplotAb)
        
        # Plot mean RoC-curves for certain experimental design with Variability-strata
        # meanROCplotV <- ggplot(data=meanROCV2, aes(x=FPR, y=meanTPR, group=strata)) + 
        #  geom_ribbon(aes(ymin=(min), ymax=(max), fill=strataClass[strata]), alpha = 0.2) + 
        #  geom_line(aes(color=strataClass[strata])) + #theme_minimal() + 
        #  labs(title=sprintf("Mean ROC-curve for %s with effect %g", plotName,q), 
        #       subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
        #       x = "False Positive Rate", y = "True Positive Rate", color="Variability strata", fill="Variability strata")+
        #  ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
        #  scale_fill_viridis_d(begin = 0.2, end = 0.6) +
        #  scale_colour_viridis_d(begin = 0.2, end = 0.6)
        #print(meanROCplotV)
        #if(savePlot == TRUE){
        #  path_save <-  sprintf("../../Result/%s/IntermediatePlots/meanROC_%s_VarStrata.pdf", saveName, saveExpDesign)
        #  ggsave(filename = path_save, plot = meanROCplotV, height = 5, width = 6)
        #  dev.off()
        #  print(meanROCplotV)
        #  rm(path_save)
        #}
        #rm(meanROCplotV)
        
        # Save strata-results for this experimental design
        for (k in 1:numberOfStrata) {
          medianAUCfinalAb <-rbind(medianAUCfinalAb,data.frame(t(colMedians(as.matrix(AUCAbundance[AUCAbundance$strata==k,])[,1:5],na.rm = T)), t(colSums(!is.na(AUCAbundance[AUCAbundance$strata==k,1:5]))), k, d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
          medianAUCfinalV <- rbind(medianAUCfinalV,data.frame(t(colMedians(as.matrix(AUCVariability[AUCVariability$strata==k,])[,1:5],na.rm = T)), t(colSums(!is.na(AUCVariability[AUCVariability$strata==k,1:5]))),k, d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
        }
        meanROCfinalAb<-rbind(meanROCfinalAb,data.frame(meanROCAb2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
        meanROCfinalV<-rbind(meanROCfinalV,data.frame(meanROCV2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
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
    dS=extraSeqSaveName[i]
    cat(sprintf("================================== d=%s =====================================\n", dD))
      
    # Creating empty initial result-matrices
    AUC = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),rep = numeric(0))
    ROC = data.frame(TPR=numeric(0), FPR=numeric(0),rep=numeric(0))
    meanROC = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0))
    genesFDR = data.frame()
    if (runStrata==T){
      AUCAbundance = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
      AUCVariability = data.frame(AUC1=numeric(0), AUC5 = numeric(0), AUCtot = numeric(0), TPR1=numeric(0), TPR5=numeric(0),run = numeric(0),strata = numeric(0))
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
        source("Resample_datasets.R", local = TRUE)
      }
      
      # Run the code for analysing DAGs
      source("Analysis_of_DAGs.R", local = TRUE)
        
      # Save the results 
      AUC[run,] <- AUCs
      ROC <- rbind(ROC,ROCs)
      meanROC<-rbind(meanROC, meanROCs)
      genesFDR <- rbind(genesFDR,genesFDRs)
      if (runStrata==T){
        AUCAbundance<-rbind(AUCAbundance, AUCsAbundance)
        AUCVariability<-rbind(AUCVariability, AUCsVariability)
        ROCAbundance <- rbind(ROCAbundance,ROCsAbundance)
        ROCVariability <- rbind(ROCVariability,ROCsVariability)
        meanROCAbundance <- rbind(meanROCAbundance,meanROCsAbundance)
        meanROCVariability <- rbind(meanROCVariability,meanROCsVariability)
        rm(ROCsAbundance, ROCsVariability, meanROCsAbundance, meanROCsVariability, AUCsAbundance,AUCsVariability)
      }
      
      rm(AUCs, ROCs, meanROCs, genesFDRs)
    }
    
    # plot individual ROCs
    # ROC$run<-as.factor(ROC$run)
    # individual_ROC_plot(ROC)
    
    colnames(meanROC)[3]<-"meanTPR"
    meanROC2<-ddply(meanROC, "FPR", summarise,
                    N    = length(meanTPR),
                    mean = mean(meanTPR),
                    min  = min(meanTPR),
                    max  = max(meanTPR))
    colnames(meanROC2)[3]<-"meanTPR"
    
    # Save results for this experimental design
    medianAUCfinal<-rbind(medianAUCfinal,data.frame(t(colMedians(as.matrix(AUC)[,1:5],na.rm = T)), t(colSums(!is.na(AUC[1:5]))), d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
    meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
    medianGenesFDR<-rbind(medianGenesFDR,data.frame(t(colMedians(as.matrix(genesFDR),na.rm = T)), t(colSums(!is.na(genesFDR))), sprintf("m=%d, d=%s",m,dD)))
    
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
      for (k in 1:numberOfStrata) {
        medianAUCfinalAb <-rbind(medianAUCfinalAb,data.frame(t(colMedians(as.matrix(AUCAbundance[AUCAbundance$strata==k,])[,1:5],na.rm = T)), t(colSums(!is.na(AUCAbundance[AUCAbundance$strata==k,1:5]))),k,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
        medianAUCfinalV <- rbind(medianAUCfinalV,data.frame(t(colMedians(as.matrix(AUCVariability[AUCVariability$strata==k,])[,1:5],na.rm = T)), t(colSums(!is.na(AUCVariability[AUCVariability$strata==k,1:5]))),k, d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
      }
      meanROCfinalAb <- rbind(meanROCfinalAb,data.frame(meanROCAb2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
      meanROCfinalV <- rbind(meanROCfinalV,data.frame(meanROCV2,d,m,m*d,sprintf("m=%d, d=%s",m,dD)))
      
      rm(ROCAbundance, ROCVariability, meanROCAbundance, meanROCVariability, meanROCAb2, meanROCV2, AUCAbundance,AUCVariability)
    }

    rm(ROC, AUC,  meanROC2, run, meanROC, genesFDR)
  }}
  
  #===================================================================================================================================
  ### Summarising results:

  colnames(medianAUCfinal)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "AUC1 Non-NA", "AUC5 Non-NA", "AUCtot Non-NA", "TPR1 Non-NA", "TPR5 Non-NA", "d", "m" ,"md","plotMD")
  colnames(meanROCfinal)<-c("FPR", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
  colnames(medianGenesFDR)<-c("Median TP count",  "Median FP count",  "Median true FDR", "TP Non-NA", "FP Non-NA", "FDR Non-NA","plotMD")
  
  # If the amount of observations is "limitNA" or less (or NA), set values to NA
  medianGenesFDR[(medianGenesFDR[,6]<=limitNA|is.na(medianGenesFDR[,6])),3]<-NA
  medianGenesFDR[(medianGenesFDR[,5]<=limitNA|is.na(medianGenesFDR[,5])),2]<-NA
  medianGenesFDR[(medianGenesFDR[,4]<=limitNA|is.na(medianGenesFDR[,4])),1]<-NA
  medianGenesFDR<-medianGenesFDR[,-c(4,5,6)]
  
  medianAUCfinal[(medianAUCfinal[,6]<=limitNA|is.na(medianAUCfinal[,6])),1]<-NA 
  medianAUCfinal[(medianAUCfinal[,7]<=limitNA|is.na(medianAUCfinal[,7])),2]<-NA 
  medianAUCfinal[(medianAUCfinal[,8]<=limitNA|is.na(medianAUCfinal[,8])),3]<-NA 
  medianAUCfinal[(medianAUCfinal[,9]<=limitNA|is.na(medianAUCfinal[,9])),4]<-NA 
  medianAUCfinal[(medianAUCfinal[,10]<=limitNA|is.na(medianAUCfinal[,10])),5]<-NA 
  medianAUCfinal<-medianAUCfinal[,-c(6:10)]
  
  #medianAUCfinal$md[medianAUCfinal$md=="5e+06"]<-"bold"
  #medianAUCfinal$md[medianAUCfinal$md==boldvalue2]<-"bold"
  #medianAUCfinal$md[medianAUCfinal$md!="bold"]<-"plain"
  
  # Save tables:
  write.csv(medianAUCfinal, file=sprintf("../../Result/%s_%s/AUC_10q%d.csv", saveName, analysis, 10*q))
  write.csv(medianGenesFDR, file=sprintf("../../Result/%s_%s/GenesFDR_10q%d.csv", saveName,analysis, 10*q))
  
  
  # # Read tables:
  # saveName="Gut2"
  # plotName="Human Gut II"
  # analysis="DESeq"
  # q=3
  # extraDesigns=T
  # extraL=3
  # savePlot=T
  # medianAUCfinal <- read.csv(file=sprintf("../../Result/%s_DESeq/AUC_10q%d.csv", saveName,10*q))[,-1]
  # medianGenesFDR <- read.csv(file=sprintf("../../Result/%s_DESeq/GenesFDR_10q%d.csv", saveName,10*q))[,-1]
  
  if (extraDesigns==T){
    HeatmapData<-data.frame(head(medianAUCfinal,-extraL), head(medianGenesFDR[1:3],-extraL)) 
  } else if (extraDesigns==F){
    HeatmapData<-data.frame(medianAUCfinal, medianGenesFDR[1:3]) 
  }
  HeatmapData$d<-as.factor(HeatmapData$d)
  HeatmapData$m<-as.factor(HeatmapData$m)
  HeatmapData$md<-as.factor(HeatmapData$md)
  
  meanROCfinal$d<-as.factor(meanROCfinal$d)
  meanROCfinal$m<-as.factor(meanROCfinal$m)
  meanROCfinal$plotMD = factor(meanROCfinal$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)

  ## for strata
  if (runStrata==T){
    medianAUCfinalAb[(medianAUCfinalAb[,6]<=limitNA|is.na(medianAUCfinalAb[,7])),1]<-NA 
    medianAUCfinalAb[(medianAUCfinalAb[,7]<=limitNA|is.na(medianAUCfinalAb[,8])),2]<-NA 
    medianAUCfinalAb[(medianAUCfinalAb[,8]<=limitNA|is.na(medianAUCfinalAb[,9])),3]<-NA 
    medianAUCfinalAb[(medianAUCfinalAb[,9]<=limitNA|is.na(medianAUCfinalAb[,10])),4]<-NA
    medianAUCfinalAb[(medianAUCfinalAb[,10]<=limitNA|is.na(medianAUCfinalAb[,11])),5]<-NA 
    medianAUCfinalAb<-medianAUCfinalAb[,-c(6:10)]
    
    medianAUCfinalV[(medianAUCfinalV[,6]<=limitNA|is.na(medianAUCfinalV[,7])),1]<-NA 
    medianAUCfinalV[(medianAUCfinalV[,7]<=limitNA|is.na(medianAUCfinalV[,8])),2]<-NA 
    medianAUCfinalV[(medianAUCfinalV[,8]<=limitNA|is.na(medianAUCfinalV[,9])),3]<-NA 
    medianAUCfinalV[(medianAUCfinalV[,9]<=limitNA|is.na(medianAUCfinalV[,10])),4]<-NA 
    medianAUCfinalV[(medianAUCfinalV[,10]<=limitNA|is.na(medianAUCfinalV[,10])),5]<-NA 
    medianAUCfinalV<-medianAUCfinalV[,-c(6:10)]
    
    colnames(medianAUCfinalAb)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "strata", "d", "m" ,"md","plotMD")
    colnames(meanROCfinalAb)<-c("FPR", "strata", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
    colnames(medianAUCfinalV)<-c("AUC1", "AUC5", "AUCtot", "TPR1", "TPR5", "strata", "d", "m" ,"md","plotMD")
    colnames(meanROCfinalV)<-c("FPR", "strata", "N", "meanTPR", "min", "max", "d", "m" ,"md","plotMD")
    
    meanROCfinalAb$d<-as.factor(meanROCfinalAb$d)
    meanROCfinalAb$m<-as.factor(meanROCfinalAb$m)
    meanROCfinalAb$plotMD = factor(meanROCfinalAb$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)
    
    meanROCfinalV$d<-as.factor(meanROCfinalV$d)
    meanROCfinalV$m<-as.factor(meanROCfinalV$m)
    meanROCfinalV$plotMD = factor(meanROCfinalV$plotMD, levels=unique(meanROCfinal$plotMD[order(meanROCfinal$m,meanROCfinal$d)]), ordered=TRUE)
  
    for (i in relations) {
      medianAUCrelationAbundance<-medianAUCfinalAb[medianAUCfinalAb$md==i,]
      medianAUCrelationAbundance<-medianAUCrelationAbundance[order(medianAUCrelationAbundance$plotMD, decreasing = TRUE),]
      medianAUCrelationAbundance<-medianAUCrelationAbundance[order(medianAUCrelationAbundance$strata),]
      medianAUCrelationAbundance<-medianAUCrelationAbundance[,-c(4,5,7,8,9)]
      write.csv(medianAUCrelationAbundance, file=sprintf("../../Result/%s_%s/AUC_Abundance_%d_10q%d.csv", saveName, analysis ,i,10*q))
      
      medianAUCrelationVariability<-medianAUCfinalV[medianAUCfinalV$md==i,]
      medianAUCrelationVariability<-medianAUCrelationVariability[order(medianAUCrelationVariability$plotMD, decreasing = TRUE),]
      medianAUCrelationVariability<-medianAUCrelationVariability[order(medianAUCrelationVariability$strata),]
      medianAUCrelationVariability<-medianAUCrelationVariability[,-c(4,5,7,8,9)]
      write.csv(medianAUCrelationVariability, file=sprintf("../../Result/%s_%s/AUC_Variability_%d_10q%d.csv", saveName, analysis, i,10*q))
    
      rm(medianAUCrelationAbundance, medianAUCrelationVariability)
    }
    
    }
  
  ### Plotting heatmaps for AUC- and TPR-values
  plot_heatmaps(HeatmapData$AUC1,"AUC", "0.01", expression(AUC[0.01]), "AUC1")
  plot_heatmaps(HeatmapData$AUC5, "AUC", "0.05", expression(AUC[0.05]), "AUC5")
  plot_heatmaps(HeatmapData$AUCtot, "AUC", "tot", expression(AUC[tot]), "AUCtot")
  plot_heatmaps(HeatmapData$TPR1, "TPR", "0.01", expression(TPR[0.01]), "TPR1")
  plot_heatmaps(HeatmapData$TPR5, "TPR", "0.05", expression(TPR[0.05]), "TPR5")
  suppressWarnings( plot_heatmaps(HeatmapData$Median.true.FDR, "True FDR", 0, "True FDR", "FDR") )
  
  ### Plot mean RoC-curves for all experimental designs
  # mean plots with set groupsize
  suppressWarnings( plot_combined_meanROCs(meanROCfinal, meanROCfinal$m, groupSize, "group size", "groupsize", meanROCfinal$d, "Sequencing depth", 0))
  # mean plots with set sequencing depth
  suppressWarnings( plot_combined_meanROCs(meanROCfinal, meanROCfinal$d, sequencingDepth, "sequencing depth", "depth", meanROCfinal$m, "Group size", 0))
  # mean plots with set groupsize and depth relation (m*d)
  suppressWarnings( plot_combined_meanROCs(meanROCfinal, meanROCfinal$md, relations, "reads in total", "relation", meanROCfinal$plotMD, "Experimental design", 0))
 
  ### Plot mean ROC-curves for strata
  if (runStrata==T){
  # mean plots with set groupsize and depth relation (m*d)
    for (strata in 1:numberOfStrata){
      class<-strataClass[strata]
      plotData=meanROCfinalAb[meanROCfinalAb$strata==strata,]
      suppressWarnings( plot_combined_meanROCs(plotData, plotData$md, relations, "reads in total", "relation", plotData$plotMD, "Experimental design", strata, sprintf("genes with %s abundance", class),"abundance"))
      plotData=meanROCfinalV[meanROCfinalV$strata==strata,]
      suppressWarnings( plot_combined_meanROCs(plotData, plotData$md, relations, "reads in total", "relation", plotData$plotMD, "Experimental design", strata, sprintf("genes with %s variability", class),"variability"))
    }
    rm(strata, plotData, class)
  }
  
  rm(group,seq, m, d)
}
rm(effect, q)