# Original datasets
Gut1 <-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt", header=T, row.names = 1)
Gut2 <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
Marine <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)

repeats=1

# For a certain number of repeats, perform the entire analysis for Gut2:
seed=1:repeats # In order to get the same results each time

for (i in 1:repeats){
  # Select dataset for resampling and downsampling
  Data=Gut2
  m=60          # Number of samples in each group (total nr samples = 2*m)
  d=2000000     # Desired sequencing depth per sample. It will not be exct
  q = 10        # Fold-change for downsampling
  f = 0.10      # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
  
  # Run the code for resampling and downsampling
  source("Resample_datasets_clean.R")
  
  # Run the code for analysing DAGs

  # ... not done
}
