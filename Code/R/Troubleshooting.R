# Troubleshooting high sequencing depths 
library(reshape2)

#=================================== Choose "good" dataset with low sequencing depth ========================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 10              # Number of samples in each group (total nr samples = 2*m)
d = 10000           # Desired sequencing depth per sample. It will not be exct
dD="10k"            # display name
q = 1.5             # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=1               # choose wich run to extract


{ # Quickly gives the case the correct name
  if (saveName == "Gut2"){
    plotName = "Human Gut II"
  } else if (saveName == "Marine"){
    plotName = "Marine"
  } else {
    sprintf("Missing name for dataset")
  }
  
  # Names for a certain dataset and name    # Results in:
  saveExpDesign = sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
  plotExpDesign = sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
  
}

# Read dataset 
DownSampledData=read.csv(sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

source("Analysis_of_DAGs_separate.R")

# Save output named after this dataset
DataLow=DownSampledData
DAGsLow=DAGs
ResDESeqLow=ResDESeq
dDLow=dD

#=================================== Choose "bad" dataset with high sequencing depth ========================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 10              # Number of samples in each group (total nr samples = 2*m)
d = 1000000         # Desired sequencing depth per sample. It will not be exct
dD="1M"             # display name
q = 1.5             # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=1               # choose wich run to extract


{ # Quickly gives the case the correct name
  if (saveName == "Gut2"){
    plotName = "Human Gut II"
  } else if (saveName == "Marine"){
    plotName = "Marine"
  } else {
    sprintf("Missing name for dataset")
  }
  
  # Names for a certain dataset and name    # Results in:
  saveExpDesign = sprintf("m%d_d%s_10q%d_f%d", m, dD, q*10, f*100)
  plotExpDesign = sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
  
}

# Read dataset 
DownSampledData=read.csv(sprintf("../../Intermediate/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

source("Analysis_of_DAGs_separate.R")

# Save output named after this dataset
DataHigh=DownSampledData
DAGsHigh=DAGs
ResDESeqHigh=ResDESeq
dDHigh=dD

# Remove non used ouptus
rm(DownSampledData, DAGs, ResDESeq, AUCs, deseqROCAUC, meanROCs, ROCs)

#========================= Dispertion and log2 fold change plots ===========================================

# check that the same DAGs have been introduced 
all(rownames(DAGsLow)==rownames(DAGsHigh))

# vector inicating if each gene is a DAG or not (DAG-indicator)
DAG=c()
for (i in 1:nrow(ResultLog2Disp)) {
  if (rownames(ResultLog2Disp[i,]) %in% rownames(DAGsLow)) {
    DAG[i]=1
  } else {DAG[i]=0}
}

# dataframe with log2 and disp for both datasets, and DAG-indicator
ResultLog2Disp=data.frame(ResDESeqLow[order(rownames(ResDESeqLow)),c(3,4)], ResDESeqHigh[order(rownames(ResDESeqHigh)),c(3,4)], DAG)
colnames(ResultLog2Disp)<- c("log2Low", "DispLow", "log2High", "DispHigh", "DAG")

# Dispertion plot
ggplot(ResultLog2Disp, aes(x=DispHigh, y=DispLow)) + 
  scale_color_viridis_d(begin = 0, end = 0.5, name=" ", labels=c("No DAG","DAG")) +
  geom_point(aes(color=factor(DAG)), size=1) +
  labs(title="Dispersion for same dataset with different seq. depths", 
       x=sprintf("Sequencing depth %s", dDHigh), y=sprintf("Sequencing depth %s", dDLow))

# log2 plot
ggplot(ResultLog2Disp, aes(x=log2High, y=log2Low)) +
  scale_color_viridis_d(begin = 0, end = 0.5, name=" ", labels=c("No DAG","DAG")) +
  geom_point(aes(color=factor(DAG)), size=1) +
  labs(title="log2 fold change for same dataset with different seq. depths", 
       x=sprintf("Sequencing depth %s", dDHigh), y=sprintf("Sequencing depth %s", dDLow))







#





















# Fulplottar 
#ResultLow=ResDESeqLow[order(rownames(ResDESeqLow)),]
#ResultHigh=ResDESeqHigh[order(rownames(ResDESeqHigh)),]
#plot(ResultHigh[,4], ResultLow[,4])
#plot(ResultHigh[,3], ResultLow[,3])

# NOT NEEDED 
#=============================== Prepare datasets for gg-plot ==============================================

# vector with rownames, ordered alphabetically/numerically
#Genes=rownames(ResDESeqLow[order(rownames(ResDESeqLow)),])  

# extract log2 for both low and high dataset
#ResultLog2=data.frame(Genes, ResDESeqLow[order(rownames(ResDESeqLow)),3], ResDESeqHigh[order(rownames(ResDESeqHigh)),3])
#colnames(ResultLog2)<- c("Gene", "Low", "High")

# extract dispertion for both low and high dataset
#ResultDisp=data.frame(Genes, ResDESeqLow[order(rownames(ResDESeqLow)),4], ResDESeqHigh[order(rownames(ResDESeqHigh)),4])
#colnames(ResultDisp)<- c("Gene", "Low", "High")

# melt both datasets and merge them 
#ResultLog2Melt=melt(ResultLog2, id="Gene", value.name = "log2", variable.name = "Type")
#ResultDispMelt=melt(ResultDisp, id="Gene", value.name = "Dispertion", variable.name = "Type")
#ResultGGFormat=merge(ResultLog2Melt,ResultDispMelt)

# remove intermediate datasets
#rm(Genes, ResultLog2, ResultDisp, ResultLog2Melt, ResultDispMelt)



