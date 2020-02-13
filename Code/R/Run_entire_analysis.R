
# Script for running all scripts, performing the entire analysis of datasets with different designs

#===================================================================================================================================
## Loading Libraries:

library(plyr)
library(viridis)
# Required packages for DAGs analysis:
library(DESeq2)
library(ggplot2)
library(pracma)

#===================================================================================================================================
## Loading data:

# Original datasets
Gut2Original <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
MarineOriginal <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)

# Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
Gut2Intermediate = Gut2Original[,colSums(Gut2Original)>=5000000]
MarineIntermediate = MarineOriginal[,colSums(MarineOriginal)>=10000000]

############## Function to Remove low counts #################################################################
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
###############################################################################################################

# Filter out genes with too low counts 
Gut2 <- remove_low_counts(Gut2Intermediate)
Marine <- remove_low_counts(MarineIntermediate)

rm(Gut2Original, Gut2Intermediate, MarineOriginal, MarineIntermediate) # remove original and intermediate datasets
#===================================================================================================================================
## General settings:

repeats = 10
seeds = 1:repeats # In order to get the same results each time

#===================================================================================================================================
## Selecting parameters and data:

Data = Marine # Gut2 or Marine
q = 2         # Fold-change for downsampling
f = 0.10      # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced

groupSize<-c(3,5,10,30,50)
sequencingDepth<-c(10000,100000,500000,1000000, 10000000) # Gut2: 5000000, Marine: 10000000

# The above sets:
# m = Number of samples in each group (total nr samples = 2*m)
# d = Desired sequencing depth per sample. It will not be exct

# Creating empty results-matrices
AUCfinal = data.frame(AUC5=numeric(0), AUC10 = numeric(0), AUCtot = numeric(0), TPR5=numeric(0), TPR10=numeric(0), sequencingDepth = character(0), groupSize = character(0))
meanROCfinal = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), sd=numeric(0),se=numeric(0),d = character(0), m = character(0))

# Looping all parameters, creating different setups
for (group in 1:length(groupSize)){
  for (seq in 1:length(sequencingDepth)) {
    m=groupSize[group]
    d=sequencingDepth[seq]
  
    { # Quickly gives the case the correct names
      if (all(dim(Data) == dim(Gut2))){
        saveName = "Gut2"
        plotName = "Human Gut II"
      } else if (all(dim(Data)==dim(Marine))){
          saveName = "Marine"
          plotName = "Marine"
      } else {
        sprintf("Missing name for dataset")
      }
      
      if (d==1e4||d==1e5||d==5e5){
        dD=d/1000
        prefix="k"
      } else if(d==1e6||d==5e6||d==10e6){
        dD=d/1000000
        prefix="M"
      } else {
        dD=d
        sprintf("wrong d")
        prefix=""
      }
      # Names for a certain dataset and name    # Results in:
      saveExpDesign = sprintf("m%d_d%d%s_q%d_f%d", m, dD, prefix, q, f*100)
      plotExpDesign = sprintf("m=%d, d=%d%s, q=%d, f=%d%%",m,dD,prefix,q,f*100)
    
      rm(dD,prefix)
      
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, saveExpDesign))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", saveExpDesign)), recursive = T)
      }
      
      if (!dir.exists(sprintf("../../Result/%s/%s", saveName, saveExpDesign))){
        dir.create(file.path("../../Result", sprintf("%s", saveName), sprintf("%s", saveExpDesign)), recursive = T)
      }
    }

    #===================================================================================================================================
    ## Run the analysis for a selected m and d:
    
    # Creating empty result-matrices
    AUC = data.frame(AUC5=numeric(0), AUC10 = numeric(0), AUCtot = numeric(0), TPR5=numeric(0), TPR10=numeric(0),seed = numeric(0))
    ROC = data.frame(TPR=numeric(0), FPR=numeric(0),seed=numeric(0))
    meanROC = data.frame(FPR=numeric(0),N=numeric(0),mean=numeric(0), sd=numeric(0),se=numeric(0))
    
    for (run in 1:repeats){
      
      selectedSeed = seeds[run]
      
      # Run the code for resampling and downsampling
      source("Resample_datasets.R")
      
      # Run the code for analysing DAGs
      savePlot = FALSE
      source("Analysis_of_DAGs.R")
    
      # Save the results 
      AUC[run,] <- AUCs
      ROC <- rbind(ROC,ROCs)
      meanROC<-rbind(meanROC, meanROCs)
      
      rm(AUCs,ROCs,meanROCs)
    }
    
    ROC$seed<-as.factor(ROC$seed)
    
    # Plot individual ROC-plots
    ROCplot <- ggplot(data=ROC, aes(x=FPR, y=TPR, group=seed)) +  geom_line(aes(color=seed)) + 
      theme(plot.title = element_text(hjust = 0.5)) +  theme_minimal() + 
      scale_color_viridis(begin = 0, end = 1, discrete=TRUE) +
      #scale_color_manual(values=colorRampPalette(brewer.pal(9, "Spectral"))(repeats))+ #c('#225EA8','#7FCDBB','#EDF8B1')) +
      labs(title=sprintf("ROC-curves for analysis of %s", plotName), 
           subtitle = sprintf("Experimental design: %s", plotExpDesign),
           colour="repeats", x = "False Positive Rate", y = "True Positive Rate") +
      xlim(0, 1)+  ylim(0, 1)
    
    print(ROCplot)
    
    savePlot=T # change when not running with 10 repeats
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/%s/individualROCs.pdf", saveName, saveExpDesign, seed)
      ggsave(filename = path_save, plot = ROCplot, height = 5, width = 6)
      dev.off()
      print(ROCplot)
    }
    
    colnames(meanROC)[3]<-"meanTPR"
    meanROC2<-ddply(meanROC, "FPR", summarise,
                    N    = length(meanTPR),
                    mean = mean(meanTPR),
                    sd   = sd(meanTPR),
                    se   = sd / sqrt(N))
    colnames(meanROC2)[3]<-"meanTPR"
    
    # Plot mean RoC-curves for certain experimental design
    meanROCplot <- ggplot(data=meanROC2, aes(x=FPR, y=meanTPR)) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd), fill="#22A88433"), alpha = 0.2) + 
      geom_line(aes(color="#22A88433")) + theme(legend.position = "none") +
      labs(title=sprintf("Mean ROC-curve for analysis of %s", plotName), 
           subtitle = sprintf("Experimental design: %s     (%s repeats)", plotExpDesign, repeats),
           x = "False Positive Rate", y = "True Positive Rate")+
      xlim(0, 1)+  ylim(0, 1) +
      scale_fill_viridis_d(begin = 0.2, end = 0.6) +
      scale_colour_viridis_d(begin = 0.2, end = 0.6)
    
    print(meanROCplot)
    
    savePlot = T # change when running with 10 repeats
    if(savePlot == TRUE){
      path_save <-  sprintf("../../Result/%s/%s/meanROC.pdf", saveName, saveExpDesign)
      ggsave(filename = path_save, plot = meanROCplot, height = 5, width = 6)
      dev.off()
      print(meanROCplot)
    }
    
    AUCfinal<-rbind(AUCfinal,c(colMeans(AUC[,1:5]),d,m))
    meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m))
    
    rm(ROCplot, meanROC, meanROCplot)
    rm(ROC, AUC,  meanROC2, run)
  }
}

#===================================================================================================================================
### Summarising results:

colnames(AUCfinal)<-c("mean AUC_{0.05}", "mean AUC_{0.10}", "mean total AUC", "mean TPR_{0.05}", "mean TPR_{0.10}", "sequencing depth", "group size" )
meanROCfinal$d<-as.factor(meanROCfinal$d)
meanROCfinal$m<-as.factor(meanROCfinal$m)
AUCfinal$`sequencing depth`<-as.factor(AUCfinal$`sequencing depth`)
AUCfinal$`group size`<-as.factor(AUCfinal$`group size`)

# Save tables:
write.csv(AUCfinal, file=sprintf("../../Result/%s/AUC.csv", saveName))

heatmap <- ggplot(AUCfinal, aes(x=AUCfinal$`group size`, y=AUCfinal$`sequencing depth`, fill=AUCfinal$`mean total AUC`)) +
                geom_tile(aes(fill = AUCfinal$`mean total AUC`)) + geom_text(aes(label = round(AUCfinal$`mean total AUC` , 2))) +
                scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
                scale_y_discrete(limits = rev(levels(as.factor(AUCfinal$`sequencing depth`)))) +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
                labs(title=sprintf("Mean total AUC-values for analysis of %s", plotName), 
                     #subtitle = sprintf("Experimental designs with group size %d     (%s repeats each)", M, repeats),
                     x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
print(heatmap)

savePlot = T # change when running with 10 repeats  
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap.pdf", saveName)
  ggsave(filename = path_save, plot = heatmap, height = 5, width = 6)
  dev.off()
  print(heatmap)}


### Plot mean RoC-curves for all experimental designs

# mean plots with set groupsize
for (group in 1:length(groupSize)){
  M=groupSize[group]
  
  meanROCplotseq <- ggplot(data=meanROCfinal[meanROCfinal$m==M,], aes(x=FPR, y=meanTPR, fill=d)) +  theme_minimal() + 
    geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd),fill = d), alpha=0.2) +
    geom_line(aes(color = d)) +
    labs(title=sprintf("Mean ROC-curves for analysis of %s", plotName), 
         subtitle = sprintf("Experimental designs with group size %d     (%s repeats each)", M, repeats),
         x = "False Positive Rate", y = "True Positive Rate",  color = "sequensing depth", fill = "sequensing depth") +
    xlim(0, 1) +  ylim(0, 1) +
    scale_fill_viridis_d(begin = 0.2, end = 0.6) +
    scale_colour_viridis_d(begin = 0.2, end = 0.6)
  
  print(meanROCplotseq)
  
  savePlot = T # change when running with 10 repeats  
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/meanROC_groupsize%d.pdf", saveName, M)
    ggsave(filename = path_save, plot = meanROCplotseq, height = 5, width = 6)
    dev.off()
    print(meanROCplotseq)}
  
  rm(meanROCplotseq, M)
}

# mean plots with set sequencing depth
for (seq in 1:length(sequencingDepth)) {
  D=sequencingDepth[seq]
  
  meanROCplotgroup <- ggplot(data=meanROCfinal[meanROCfinal$d==D,], aes(x=FPR, y=meanTPR, fill=m)) +  theme_minimal() + 
    geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd),fill = m), alpha=0.2) +
    geom_line(aes(color = m)) +
    labs(title=sprintf("Mean ROC-curves for analysis of %s", plotName), 
         subtitle = sprintf("Experimental designs with sequencing depth %d     (%s repeats each)", D, repeats),
         x = "False Positive Rate", y = "True Positive Rate",  color = "group size", fill = "group size") +
    xlim(0, 1) +  ylim(0, 1) +
    scale_fill_viridis_d(begin = 0.2, end = 0.6) +
    scale_colour_viridis_d(begin = 0.2, end = 0.6)
  
  print(meanROCplotgroup)
  
  savePlot=T # change when running with 10 repeats
  if(savePlot == TRUE){
    path_save <-  sprintf("../../Result/%s/meanROC_depth%d.pdf", saveName, D)
    ggsave(filename = path_save, plot = meanROCplotgroup, height = 5, width = 6)
    dev.off()
    print(meanROCplotgroup)}
  
  rm(meanROCplotgroup, D)
}

rm(group,seq, repeats)

