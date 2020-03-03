# Troubleshooting high sequencing depths 

#=================================== Choose "good" dataset with low sequencing depth ========================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 3               # Number of samples in each group (total nr samples = 2*m)
d = 500000          # Desired sequencing depth per sample. It will not be exct
dD="500k"           # display name
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
DownSampledData=read.csv(sprintf("../../Intermediate/%s_sameseed/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/%s_sameseed/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

source("Analysis_of_DAGs_separate.R")

# Save output named after this dataset
DataLow=DownSampledData
DAGsLow=DAGs
ResDESeqLow=ResDESeq
dDLow=dD

#=================================== Choose "bad" dataset with high sequencing depth ========================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 3               # Number of samples in each group (total nr samples = 2*m)
d = 10000000        # Desired sequencing depth per sample. It will not be exct
dD="10M"            # display name
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
DownSampledData=read.csv(sprintf("../../Intermediate/%s_sameseed/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/%s_sameseed/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

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

# dataframe with log2 and disp for both datasets
ResultLog2Disp=data.frame(ResDESeqLow[order(rownames(ResDESeqLow)),c(3,4)], ResDESeqHigh[order(rownames(ResDESeqHigh)),c(3,4)])
colnames(ResultLog2Disp)<- c("log2Low", "DispLow", "log2High", "DispHigh")

# vector inicating if each gene is a DAG or not (DAG-indicator)
DAG=c()
for (i in 1:nrow(ResultLog2Disp)) {
  if (rownames(ResultLog2Disp[i,]) %in% rownames(DAGsLow)) {
    DAG[i]=1
  } else {DAG[i]=0}
}

# add DAG-indicator to dataframe
ResultLog2Disp=data.frame(ResultLog2Disp, DAG)


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


#================================= Investigate most significant genes =====================================

##################################
# Prepare "low" dataset:
#DAG-indicator
DAG=c()
for (i in 1:nrow(ResDESeqLow)) {
  if (rownames(ResDESeqLow[i,]) %in% rownames(DAGsLow)) {
    DAG[i]=1
  } else {DAG[i]=0}
}
# Dataframe with the names of the genes, p-values and DAG-indicator (ordered with increasing p-value)
GenesLow=data.frame(row.names = rownames(ResDESeqLow), ResDESeqLow[,1] , DAG)
colnames(GenesLow)<-c("p-value","DAG")

# Prepare "high" dataset
# DAG-indicator
DAG=c()
for (i in 1:nrow(ResDESeqHigh)) {
  if (rownames(ResDESeqHigh[i,]) %in% rownames(DAGsHigh)) {
    DAG[i]=1
  } else {DAG[i]=0}
}
# Dataframe with the names of the genes, p-values and DAG-indicator (ordered with increasing p-value)
GenesHigh=data.frame(row.names = rownames(ResDESeqHigh), ResDESeqHigh[,1] , DAG)
colnames(GenesHigh)<-c("p-value","DAG")

rm(DAG, ResDESeqHigh, ResDESeqLow)
##################################


# Extract top 100 genes 
GenesLow100=head(GenesLow, 100)
GenesHigh100=head(GenesHigh, 100)

# Does "high" have less TP than "low"? 
sum(GenesHigh100[,2])<sum(GenesLow100[,2])

# vector with value 1 if TP in low is also found in high and value 0 if TP is not found in high 
TPBoth=c()
for (i in 1:sum(GenesLow100[,2])) {
  if (rownames(GenesLow100[GenesLow100[,2]==1,])[i] %in% rownames(GenesHigh100[GenesHigh100[,2]==1,])) {
    TPBoth[i]=1
  } else {TPBoth[i]=0}
}

# How many TP are found in low but not in high:
sum(TPBoth==0)

# Names of these genes
TPDiffNames=rownames(GenesLow100[GenesLow100[,2]==1,])[TPBoth==0]

#  The downsampled data with only the TP genes that differ between high and low
TPDiffHigh=DataHigh[TPDiffNames,]
TPDiffLow=DataLow[TPDiffNames,]












# vector with value 1 if FP in high is also found in low and value 0 if FP is not found in low 
FPBoth=c()
for (i in 1:(100-sum(GenesHigh100[,2]))) {
  if (rownames(GenesHigh100[GenesHigh100[,2]==0,])[i] %in% rownames(GenesLow100[GenesLow100[,2]==0,])) {
    FPBoth[i]=1
  } else {FPBoth[i]=0}
}

# How many FP are found in high but not in low:
sum(FPBoth==0)

# Names of these genes
FPDiffNames=rownames(GenesHigh100[GenesHigh100[,2]==0,])[TPBoth==0]

#  The downsampled data with only the FP genes that differ between high and low
FPDiffHigh=DataHigh[FPDiffNames,]
FPDiffLow=DataLow[FPDiffNames,]








## Sanity check:
sum(GenesHigh100[,2]==1) # TP
sum(GenesHigh100[,2]==0) # FP

sum(GenesLow100[,2]==1) # TP
sum(GenesLow100[,2]==0) # FP

# How many TP are found in low but not in high:
sum(TPBoth==0)

# How many TP are found in low and in high:
sum(TPBoth==1)

# How many FP are found in high but not in low:
sum(FPBoth==0)

# How many FP are found in high and in low:
sum(FPBoth==1)
