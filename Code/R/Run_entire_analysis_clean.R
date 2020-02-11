# Original datasets
#Gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
#Gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
#Marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)

library(plyr)

colorScale9<-c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
color1=colorScale9[7]

# CREATE ALL FOLDERS REQUIRED!


# General settings
repeats = 2
seeds = 1:repeats # In order to get the same results each time


# Experimental design 1 for selected data
saveExpDesign = "m60_d2e6_q2_f010" # or name it otherwise. NOTE: f0.10 doesn't work
plotExpDesign = "m = 60, d = 2000000, q = 2, f = 0.10"
m = 60        # Number of samples in each group (total nr samples = 2*m)
d = 2000000   # Desired sequencing depth per sample. It will not be exct
q = 2        # Fold-change for downsampling
f = 0.10      # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced


### GUT 2 ###
Data = read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1) # Gut2
saveName = "Gut2"
plotName = "Human Gut II"


# For a certain number of repeats, perform the entire analysis for a selected experimental design for Gut2:
AUC = data.frame(AUC5=numeric(0), AUC10 = numeric(0), AUCtot = numeric(0), seed = numeric(0))
ROC = data.frame(TPR=numeric(0), FPR=numeric(0),seed=numeric(0))
meanROC = data.frame(FPR=numeric(0),N=numeric(0),mean=numeric(0), sd=numeric(0),se=numeric(0))


for (run in 1:repeats){
  
  selectedSeed = seeds[run]
  
  # Run the code for resampling and downsampling
  source("Resample_datasets_clean.R")
  
  # Run the code for analysing DAGs
  savePlot = FALSE
  source("Analysis_of_DAGs_clean.R")

  # Save the results 
  AUC[run,] <- AUCs
  ROC <- rbind(ROC,ROCs)
  meanROC<-rbind(meanROC, meanROCs)
  #rm(AUCs,ROCsmeanROCs)
}

#rm(plotExpDesign,saveExpDesign)

colnames(meanROC)[3]<-"meanTPR"
meanROC2<-ddply(meanROC, "FPR", summarise,
              N    = length(meanTPR),
              mean = mean(meanTPR),
              sd   = sd(meanTPR),
              se   = sd / sqrt(N)
)


# Plot individual ROC-plots

ROCplot <- ggplot(data=ROC, aes(x=FPR, y=TPR), color=seed) +  geom_line() + 
  theme(plot.title = element_text(hjust = 0.5)) +  theme_minimal() + 
  scale_color_manual(values=c('#7FCDBB','#225EA8')) +
  labs(title=sprintf("ROC-curves for analysis of %s", plotName), 
       subtitle = sprintf("Experimental design: %s", plotExpDesign),
       x = "False Positive Rate", y = "True Positive Rate")

print(ROCplot)

if(savePlot == TRUE){
  path_save <-  sprintf("../../Result/%s/%s/ROC_seed%d.pdf", saveName, saveExpDesign, seed)
  ggsave(filename = path_save, plot = ROCplot, height = 5, width = 6)
  dev.off()
  print(ROCplot)
}

# Plot mean AUC-values/RoC-curves
#savePlot = FALSE

meanROCplot <- ggplot(data=meanROC2, aes(x=FPR, y=mean)) +  theme_minimal() + 
  geom_ribbon(aes(ymin=(mean-sd), ymax=(mean+sd)), alpha = 0.2, fill = color1) + geom_line(color=color1, size=1) +
  labs(title=sprintf("Mean ROC-curve for analysis of %s", plotName), 
       subtitle = sprintf("Experimental design: %s", plotExpDesign),
       x = "False Positive Rate", y = "True Positive Rate")

print(meanROCplot)
