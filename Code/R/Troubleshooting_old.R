# Troubleshooting high sequencing depths 

#=================================== Choose dataset to troubleshoot ===========================================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 50              # Number of samples in each group (total nr samples = 2*m)
d = 10000000        # Desired sequencing depth per sample. It will not be exct
dD="10M"            # display name
q = 1.5   #1.5, 2.5          # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=10    #7,10           # choose wich run to extract


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
DownSampledData=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

# =================================== Source Analysis of DAGs ==================================================

source("Analysis_of_DAGs_separate.R")

#======================================Plot TP and FP in significance-order ====================================

PR=c()
for (i in 1:nrow(DownSampledData)) {
  
  if (BinTP[i]==0){
    PR[i]="FP"
  } else {
    PR[i]="TP"
  }
}
x=seq(1,nrow(DownSampledData), by=1)
y=c(rep(1, nrow(DownSampledData)))
HeatmapPR=data.frame(x,y, PR)

ggplot(HeatmapPR, aes(x=x, y=y, fill=PR)) +
  geom_tile(aes(fill = PR)) +
  #scale_fill_manual(values = c("red", "blue")) + 
  scale_fill_viridis_d(begin = 0, end = 1, alpha = 0.5) +  
  labs(title="TP and FP in significance order", x = "Genes (from most to least significant)" , y = "  ") #, fill = "AUC-values") 

#======================================= Dispersion ============================================================
# Look att dispersion
summary(ResDESeq[,4])


#=================================Summary of log2 fold change for TP and FP =====================================
# check summary for TP and FP

# log2 fold change <0, TP 
ResultLTrue=ResDESeq[rownames(ResDESeq) %in% rownames(DAGs) & ResDESeq[,3]<0 & ResDESeq[,2]<0.05 ,]
summary(ResultLTrue[,3])

# log2 fold change <0, FP  
ResultLFalse=ResDESeq[!rownames(ResDESeq) %in% rownames(DAGs) & ResDESeq[,3]<0 & ResDESeq[,2]<0.05 ,]
summary(ResultLFalse[,3])

# log2 fold change >0, TP  
ResultBTrue=ResDESeq[rownames(ResDESeq) %in% rownames(DAGs) & ResDESeq[,3]>0 & ResDESeq[,2]<0.05 ,]
summary(ResultBTrue[,3])

# log2 fold change >0, FP  
ResultBFalse=ResDESeq[!rownames(ResDESeq) %in% rownames(DAGs) & ResDESeq[,3]>0 & ResDESeq[,2]<0.05 ,]
summary(ResultBFalse[,3])


# ============================================= Check for zeros ================================================= 
# Total number of zeros in the dataset
sum(DownSampledData==0)

# create vector containg the number of zeroes per gene
zerosPerGene=c()
for (i in 1:nrow(DownSampledData)) {
  zerosPerGene[i]=sum(DownSampledData[i,]==0)
}

summary(zerosPerGene)

# See if the gene with maximum zeros is among the introduced DAGs (if yes, returns at what index)
match(row.names(DownSampledData[match(92, zerosPerGene),]), rownames(DAGs))


