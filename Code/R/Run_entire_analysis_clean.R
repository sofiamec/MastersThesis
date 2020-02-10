# Original datasets
#Gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
#Gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
#Marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)

library(plyr)


# Function for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

# CREATE ALL FOLDERS REQUIRED!



# General settings
repeats = 3
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

  # Save the results (currently for edgeR)
  AUC[run,] <- AUCs
  ROC <- rbind(ROC,ROCs)
  meanROC<-rbind(meanROC, meanROCs)
  #rm(AUCs,ROCsmeanROCs)
}

#rm(plotExpDesign,saveExpDesign)

meanROC2<-ddply(meanROC, "FPR", summarise,
              N    = length(mean),
              mean = mean(mean),
              sd   = sd(mean),
              se   = sd / sqrt(N)
)

# Plot mean AUC-values/RoC-curves
#savePlot = FALSE

meanROCplot <- ggplot(data=meanROC2, aes(x=FPR, y=mean)) +  geom_line() +  theme_minimal() + 
  labs(title=sprintf("Mean ROC-curve for analysis of %s", plotName), 
       subtitle = sprintf("Experimental design: %s", plotExpDesign),
       x = "False Positive Rate", y = "True Positive Rate")

print(meanROCplot)
