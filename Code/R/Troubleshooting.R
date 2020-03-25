# Troubleshooting high sequencing depths
library(xtable)

analysisDESeq2=F        # Set to T if the analysis should be done with DESeq2, F is it should be done with OGLM

#=================================== Choose "good" dataset with low sequencing depth ========================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 10              # Number of samples in each group (total nr samples = 2*m)
d = 500000          # Desired sequencing depth per sample. It will not be exct
dD="500k"           # display name
q = 1.5             # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=10              # choose wich run to extract


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
ResDAGsAnalysisLow=ResDAGsAnalysis
dDLow=dD

#=================================== Choose "bad" dataset with high sequencing depth ========================
#saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
#m = 10              # Number of samples in each group (total nr samples = 2*m)
d = 10000000         # Desired sequencing depth per sample. It will not be exct
dD="10M"             # display name
#q = 1.5             # Fold-change for downsampling
#f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
#run=                # choose wich run to extract


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
ResDAGsAnalysisHigh=ResDAGsAnalysis
dDHigh=dD

# Remove non used ouptus
rm(DownSampledData,DAGs,ResDAGsAnalysis,AUCs,deseqROCAUC,meanROCs,ROCs,BinTP,dD,i,matchDESeq,d,plotExpDesign,saveExpDesign)


#========================= Dispertion and log2 fold change plots ===========================================

# check that the same DAGs have been introduced 
all(rownames(DAGsLow)==rownames(DAGsHigh))

if(analysisDESeq2==T){
  # dataframe with log2 and disp for both datasets
  ResultLog2Disp=data.frame( ResDAGsAnalysisLow[order(rownames( ResDAGsAnalysisLow)),c(3,4)],  ResDAGsAnalysisHigh[order(rownames( ResDAGsAnalysisHigh)),c(3,4)])
  colnames(ResultLog2Disp)<- c("log2Low", "DispLow", "log2High", "DispHigh")
  
  # vector indicating if each gene is a DAG or not (DAG-indicator)
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
}

#================================= Investigate most significant genes =====================================

#####################################################################
# Prepare "low" dataset:
#DAG-indicator
DAG=c()
for (i in 1:nrow( ResDAGsAnalysisLow)) {
  if (rownames( ResDAGsAnalysisLow[i,]) %in% rownames(DAGsLow)) {
    DAG[i]=1
  } else {DAG[i]=0}
}
# Dataframe with the names of the genes, p-values and DAG-indicator (ordered with increasing p-value)
GenesLow=data.frame(row.names = rownames( ResDAGsAnalysisLow),  ResDAGsAnalysisLow[,1] , DAG)
colnames(GenesLow)<-c("p-value","DAG")

# Prepare "high" dataset
# DAG-indicator
DAG=c()
for (i in 1:nrow( ResDAGsAnalysisHigh)) {
  if (rownames( ResDAGsAnalysisHigh[i,]) %in% rownames(DAGsHigh)) {
    DAG[i]=1
  } else {DAG[i]=0}
}
# Dataframe with the names of the genes, p-values and DAG-indicator (ordered with increasing p-value)
GenesHigh=data.frame(row.names = rownames( ResDAGsAnalysisHigh),  ResDAGsAnalysisHigh[,1] , DAG)
colnames(GenesHigh)<-c("p-value","DAG")

rm(DAG)
#####################################################################


# Extract top 100 genes 
GenesLow100=head(GenesLow, 100)
GenesHigh100=head(GenesHigh, 100)

# Does "high" have less TP than "low"? 
sum(GenesHigh100[,2])<sum(GenesLow100[,2])


########################### TP #######################################
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

#  The downsampled data with only the TP genes that are in low but not in high
TPDiffHigh=DataHigh[TPDiffNames,]
TPDiffLow=DataLow[TPDiffNames,]

# The position of these genes in the ordered result for the high dataset
TPPosHigh=c()
for (i in 1:length(TPDiffNames)) {
  TPPosHigh[i]=which(rownames( ResDAGsAnalysisHigh)==TPDiffNames[i])
}

# Dataframes with row means for the different groups and dispersion, adjusted p-values and position for the 
# TP-genes that appear in low but not in high
TPDiffMeansHigh=data.frame(row.names = TPDiffNames, Mean1=as.integer(rowMeans(TPDiffHigh[,c(1:m)])), 
                           Mean2=as.integer(rowMeans(TPDiffHigh[,c((m+1):(2*m))])),
                           Dispersion= ResDAGsAnalysisHigh[TPDiffNames,4], adjusted.p.value= ResDAGsAnalysisHigh[TPDiffNames,2],
                           position=TPPosHigh)

TPDiffMeansLow=data.frame(row.names = TPDiffNames, Mean1=as.integer(rowMeans(TPDiffLow[,c(1:m)])), 
                          Mean2=as.integer(rowMeans(TPDiffLow[,c((m+1):(2*m))])),
                          Dispersion= ResDAGsAnalysisLow[TPDiffNames,4], adjusted.p.value= ResDAGsAnalysisLow[TPDiffNames,2])


########################### FP #######################################
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
FPDiffNames=rownames(GenesHigh100[GenesHigh100[,2]==0,])[FPBoth==0]

#  The downsampled data with only the FP genes that are in high but not in low
FPDiffHigh=DataHigh[FPDiffNames,]
FPDiffLow=DataLow[FPDiffNames,]

# The position of these genes in the ordered result for the high dataset
FPPosHigh=c()
for (i in 1:length(FPDiffNames)) {
  FPPosHigh[i]=which(rownames( ResDAGsAnalysisHigh)==FPDiffNames[i])
}

# Dataframes with row means for the different groups and dispersion, adjusted p-values and position for the 
# FP-genes that appear in high but not in low
FPDiffMeansHigh=data.frame(row.names = FPDiffNames, Mean1=as.integer(rowMeans(FPDiffHigh[,c(1:m)])), 
                           Mean2=as.integer(rowMeans(FPDiffHigh[,c((m+1):(2*m))])),
                           Dispersion= ResDAGsAnalysisHigh[FPDiffNames,4], adjusted.p.value= ResDAGsAnalysisHigh[FPDiffNames,2],
                           Position=FPPosHigh)

FPDiffMeansLow=data.frame(row.names = FPDiffNames, Mean1=as.integer(rowMeans(FPDiffLow[,c(1:m)])), 
                          Mean2=as.integer(rowMeans(FPDiffLow[,c((m+1):(2*m))])),
                          Dispersion= ResDAGsAnalysisLow[FPDiffNames,4], adjusted.p.value= ResDAGsAnalysisLow[FPDiffNames,2])

#print(xtable(cbind(FPDiffMeansHigh, FPDiffMeansLow)))
#print(xtable(cbind(TPDiffMeansHigh, TPDiffMeansLow)))




## Sanity check:
#sum(GenesHigh100[,2]==1) # TP
#sum(GenesHigh100[,2]==0) # FP

#sum(GenesLow100[,2]==1) # TP
#sum(GenesLow100[,2]==0) # FP

# How many TP are found in low but not in high:
#sum(TPBoth==0)

# How many TP are found in low and in high:
#sum(TPBoth==1)

# How many FP are found in high but not in low:
#sum(FPBoth==0)

# How many FP are found in high and in low:
#sum(FPBoth==1)