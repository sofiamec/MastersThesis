
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

#===================================================================================================================================
## Function to Remove low counts 
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

## Function for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

#===================================================================================================================================
# Filter out genes with too low counts 
Gut2 <- remove_low_counts(Gut2Intermediate)
Marine <- remove_low_counts(MarineIntermediate)

rm(Gut2Original, Gut2Intermediate, MarineOriginal, MarineIntermediate) # remove original and intermediate datasets

####################################################################################################################################
#===================================================================================================================================
#                      Only change these parameters for different results! 
#                             (the rest of the code should adjust)
#===================================================================================================================================
## Selecting parameters and data:

repeats = 10
savePlot = T

Data = Gut2                                                 # Gut2 or Marine
effectsizes=c(1.5,1.8,2,2.5,4)                              # q 
# remove q from seeds when value is fixed!
groupSize<-c(3,5,10,30,50)                                  # m
sequencingDepth<-c(10000,100000,500000,1000000,5000000)     # d,  Gut2: 5000000, Marine: 5000000 and 10000000
sequencingDepthName<-c("10k","100k","500k","1M", "5M")      # dD, Gut2: 5M, Marine: 10M
boldvalue2="0"                                                # Gut2: "0", Marine: "50000000"

# The above sets:
# q = Fold-change for downsampling
# m = Number of samples in each group (total nr samples = 2*m)
# d = Desired sequencing depth per sample
# boldvalue2 = What experimental design, apart from "5M", that should be bold in heatmaps   

#===================================================================================================================================
#===================================================================================================================================
####################################################################################################################################

## old:
#seeds = 1:(repeats*length(groupSize)*length(sequencingDepth)*length(effectsizes)) # In order to get the same results each time
#selectedSeed = seeds[run] # used to be in inner for-llop

set.seed(100)

for (effect in 1:length(effectsizes)) { # looping over q
  q=effectsizes[effect]

f = 0.10      # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced

# Creating empty results-matrices
meanAUCfinal = data.frame(AUC5=numeric(0), AUC10 = numeric(0), AUCtot = numeric(0), TPR5=numeric(0), TPR10=numeric(0), d = character(0), m = character(0))
meanROCfinal = data.frame(FPR=numeric(0),N=numeric(0),meanTPR=numeric(0), sd=numeric(0),se=numeric(0),d = character(0), m = character(0))

### Looping all parameters, creating different setups
for (group in 1:length(groupSize)){ # looping over m
  m=groupSize[group]
  cat(sprintf("================================== m=%d =======================================\n", m))
  for (seq in 1:length(sequencingDepth)) { # looping over d
    d=sequencingDepth[seq]
    cat(sprintf("================================= d=%d ====================================\n", d))
  
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
      saveExpDesign = sprintf("m%d_d%d%s_10q%d_f%d", m, dD, prefix, q*10, f*100)
      plotExpDesign = sprintf("m=%d, d=%d%s, q=%g, f=%d%%",m,dD,prefix,q,f*100)
    
      rm(dD,prefix)
      
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, saveExpDesign))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", saveExpDesign)), recursive = T)
      }
      
      if (!dir.exists(sprintf("../../Result/%s/IntermediatePlots", saveName))){
        dir.create(file.path("../../Result", sprintf("%s", saveName), "/IntermediatePlots"), recursive = T)
      }
    }

    #===================================================================================================================================
    ## Run the analysis for a selected m and d:
    
    # Creating empty result-matrices
    AUC = data.frame(AUC5=numeric(0), AUC10 = numeric(0), AUCtot = numeric(0), TPR5=numeric(0), TPR10=numeric(0),rep = numeric(0))
    ROC = data.frame(TPR=numeric(0), FPR=numeric(0),rep=numeric(0))
    meanROC = data.frame(FPR=numeric(0),N=numeric(0),mean=numeric(0), sd=numeric(0),se=numeric(0))
    
    for (run in 1:repeats){
      cat(sprintf("Repeat %d\n", run))
      
      # Run the code for resampling and downsampling
      source("Resample_datasets.R")
      
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
      labs(title=sprintf("ROC-curves for %s", plotName), 
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
                    sd   = sd(meanTPR),
                    se   = sd / sqrt(N))
    colnames(meanROC2)[3]<-"meanTPR"
    
    # Plot mean RoC-curves for certain experimental design
    meanROCplot <- ggplot(data=meanROC2, aes(x=FPR, y=meanTPR)) +  theme_minimal() + 
      geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd), fill="#22A88433"), alpha = 0.2) + 
      geom_line(aes(color="#22A88433")) + theme(legend.position = "none") +
      labs(title=sprintf("Mean ROC-curve for %s", plotName), 
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
    
    meanAUCfinal<-rbind(meanAUCfinal,c(colMeans(AUC[,1:5]),d,m))
    meanROCfinal<-rbind(meanROCfinal,data.frame(meanROC2,d,m))
    
    rm(ROCplot, meanROC, meanROCplot)
    rm(ROC, AUC,  meanROC2, run)
  }
}

#===================================================================================================================================
### Summarising results:

meanAUCfinal<-data.frame(meanAUCfinal,(meanAUCfinal[,6]*meanAUCfinal[,7]))
colnames(meanAUCfinal)<-c("AUC5", "AUC10", "AUCtot", "TPR5", "TPR10", "d", "m" ,"md")

meanAUCfinal$md[meanAUCfinal$md=="5000000"]<-"bold"
meanAUCfinal$md[meanAUCfinal$md==boldvalue2]<-"bold"
meanAUCfinal$md[meanAUCfinal$md!="bold"]<-"plain"

meanAUCfinal$d<-as.factor(meanAUCfinal$d)
meanAUCfinal$m<-as.factor(meanAUCfinal$m)
meanAUCfinal$md<-as.factor(meanAUCfinal$md)
meanROCfinal$d<-as.factor(meanROCfinal$d)
meanROCfinal$m<-as.factor(meanROCfinal$m)

# Save tables:
write.csv(meanAUCfinal, file=sprintf("../../Result/%s/AUC_10q%d.csv", saveName,10*q))



# heatmaps for AUC and TPR at FPR 0.5, 0.10 and 1
heatmapAUCtot <- ggplot(meanAUCfinal, aes(x=meanAUCfinal$m, y=meanAUCfinal$d, fill=AUCtot)) +
  geom_tile(aes(fill = AUCtot)) + geom_text(aes(label = round2(AUCtot, 2), fontface=md)) +
  scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
  scale_y_discrete(limits = rev(levels(as.factor(meanAUCfinal$d)))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
  labs(title=sprintf("Mean total AUC-values for %s", plotName), 
       x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
print(heatmapAUCtot)
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap_AUCtot_10q%d.pdf", saveName,q*10)
  ggsave(filename = path_save, plot = heatmapAUCtot, height = 5, width = 6)
  dev.off()
  print(heatmapAUCtot)}

heatmapAUC5 <- ggplot(meanAUCfinal, aes(x=meanAUCfinal$m, y=meanAUCfinal$d, fill=AUC5)) +
  geom_tile(aes(fill = AUC5)) + geom_text(aes(label = round2(AUC5, 2), fontface=md)) +
  scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
  scale_y_discrete(limits = rev(levels(as.factor(meanAUCfinal$d)))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
  labs(title=sprintf("Mean AUC-values at FPR 0.05 for %s", plotName), 
       x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
print(heatmapAUC5)
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap_AUC5_10q%d.pdf", saveName,10*q)
  ggsave(filename = path_save, plot = heatmapAUC5, height = 5, width = 6)
  dev.off()
  print(heatmapAUC5)}

heatmapAUC10 <- ggplot(meanAUCfinal, aes(x=meanAUCfinal$m, y=meanAUCfinal$d, fill=AUC10)) +
  geom_tile(aes(fill = AUC10)) + geom_text(aes(label = round2(AUC10, 2), fontface=md)) +
  scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
  scale_y_discrete(limits = rev(levels(as.factor(meanAUCfinal$d)))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
  labs(title=sprintf("Mean AUC-values at FPR 0.10 for %s", plotName), 
       x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "AUC-values") 
print(heatmapAUC10)
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap_AUC10_10q%d.pdf", saveName,q*10)
  ggsave(filename = path_save, plot = heatmapAUC10, height = 5, width = 6)
  dev.off()
  print(heatmapAUC10)}

heatmapTPR5 <- ggplot(meanAUCfinal, aes(x=meanAUCfinal$m, y=meanAUCfinal$d, fill=TPR5)) +
  geom_tile(aes(fill = TPR5)) + geom_text(aes(label = round2(TPR5, 2), fontface=md)) +
  scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
  scale_y_discrete(limits = rev(levels(as.factor(meanAUCfinal$d)))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
  labs(title=sprintf("Mean TPR-values at FPR 0.05 for %s", plotName), 
       x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "TPR-values") 
print(heatmapTPR5)
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap_TPR5_10q%d.pdf", saveName,q*10)
  ggsave(filename = path_save, plot = heatmapTPR5, height = 5, width = 6)
  dev.off()
  print(heatmapTPR5)}

heatmapTPR10 <- ggplot(meanAUCfinal, aes(x=meanAUCfinal$m, y=meanAUCfinal$d, fill=TPR10)) +
  geom_tile(aes(fill = TPR10)) + geom_text(aes(label = round2(TPR10, 2), fontface=md)) +
  scale_fill_viridis_c(begin = 0, end = 1, alpha = 0.5) +  
  scale_y_discrete(limits = rev(levels(as.factor(meanAUCfinal$d)))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background=element_rect(fill = "white") ) +
  labs(title=sprintf("Mean TPR-values at FPR 0.10 for %s", plotName), 
       x = "Group size", y = "Sequencing depth",  color = "sequensing depth", fill = "TPR-values") 
print(heatmapTPR10)
if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/heatmap_TPR10_10q%d.pdf", saveName,q*10)
  ggsave(filename = path_save, plot = heatmapTPR10, height = 5, width = 6)
  dev.off()
  print(heatmapTPR10)}


### Plot mean RoC-curves for all experimental designs

# mean plots with set groupsize
for (group in 1:length(groupSize)){
  M=groupSize[group]
  
  meanROCplotgroup <- ggplot(data=meanROCfinal[meanROCfinal$m==M,], aes(x=FPR, y=meanTPR, fill=d)) +  theme_minimal() + 
    geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd),fill = d), alpha=0.2) +
    geom_line(aes(color = d)) +
    labs(title=sprintf("Mean ROC-curves for %s", plotName), 
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
    geom_ribbon(aes(ymin=(meanTPR-sd), ymax=(meanTPR+sd),fill = m), alpha=0.2) +
    geom_line(aes(color = m)) +
    labs(title=sprintf("Mean ROC-curves for %s", plotName), 
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

rm(group,seq)

}

rm(repeats)
