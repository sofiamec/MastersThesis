# Troubleshooting high sequencing depths 


#=================================== Choose dataset to troubleshoot ===========================================
saveName = "Marine" # Choose dataset. Ex: "Gut2" or "Marine"
m = 50              # Number of samples in each group (total nr samples = 2*m)
d = 10000000        # Desired sequencing depth per sample. It will not be exct
q = 1.5             # Fold-change for downsampling
f = 0.10            # Desired total fraction of genes to be downsampled. It will not be exact. The effects will be balanced
run=7               # choose wich run to extract

{ # Quickly gives the case the correct name
  if (saveName == "Gut2"){
    plotName = "Human Gut II"
  } else if (saveName == "Marine"){
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
}

# Read dataset 
DownSampledData=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DownSampledData_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)
DAGs=read.csv(sprintf("../../Intermediate/test_qvalues/%s/%s/DAGs_run%d.csv", saveName, saveExpDesign, run), header = T, row.names = 1)

# ================================================================================================================================

# run dataset in Analysis of DAGs
source("Analysis_of_DAGs_separate.R")

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


