# Exploration of data

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(xtable)
library(tidyverse)

#======================================== Functions ======================================#

## Reshaping dataset
# For a given dataset, this function reshapes the data into a format required by gg-plot
# The dataset has to have unique rownames corresponding to the different genes
# Input: Data = the original dataset to be reshaped (rownames required)
# Output: DataGG = the reshaped dataset with 3 columns (Gene, SampleID and GeneCount)
reshaping_dataset=function(Data){
  names(Data)[1]<-"Gene"
  DataGG<-melt(Data,id=c("Gene"))
  colnames(DataGG)<-c("Gene","SampleID","GeneCount")
  return(DataGG)
}

## Compute low counts
# For a given dataset, this function computes genes which should be removed when filtering.
# It also prints how many genes will be removed respectively kept.
# input: Data = the data to be filtered
# output: r = a logical vector with 1 if a gene should be removed, and 0 if it should be kept. 
compute_low_counts = function(Data){
  
  a=rowSums(Data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(Data)) {
    b[i] <- sum(Data[i,]==0)/ncol(Data)>=0.75
    r[i] <- a[i]+b[i]
  }
  cat(sprintf("Number of genes with low counts:        %s\n", sum(r!=0)))
  cat(sprintf("Number of genes with acceptable counts: %s\n", sum(r==0)))
  return(r)
}

# gene_boxplot
# Displays and saves a figure of how genes are distributed in all samples of the given dataset
# Input: DataGG = the reshaped dataset to be plotted
#        title = a string of what should be displayed as the figure title. Ex: "Number of reads in Human Gut I"
#        colorScale = a discrete scale of colours used for filling the boxplots. Ex: ColorScale150
#        savePlot = TRUE or FALSE, statement indicates if the plot should be saved or not
#        saveName = character with name of the dataset as it will be saved. Ex: "Gut1"
gene_boxplot = function(DataGG, title, colorScale, savePlot, saveName){
  p1 <- ggplot(DataGG,aes(x=SampleID,y=GeneCount)) + 
    geom_boxplot(aes(fill=SampleID)) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),legend.position = "none") + 
    labs(title=title) + 
    scale_fill_manual(values=colorScale) 
  
  print(p1) 
  # Saving plots:
  if(savePlot == TRUE){
    path_save1 <- str_c("../../Result/", saveName,"/", saveName, "_gene_boxplot.pdf")
    ggsave(filename = path_save1, plot = p1, height = 5, width = 12)
    dev.off()
    print(p1) 
  }
}

# sequencing_depth_histogram
# Displays and saves a figure of how the sequencing depth of samples are distributed in a dataset
# Input: DataGG = the reshaped dataset to be plotted
#        bins = number of partitions for the histogram
#        title = a string of what should be displayed as the figure title. Ex: "Sequencing Depth in Human Gut I"
#        fillcolor = colour code for filling the histogram bars. Ex: colour1
#        bordercolor = colour code for filling the bar borders. Ex: "black"
#        savePlot = TRUE or FALSE, statement indicates if the plot should be saved or not
#        saveName = character with name of the dataset as it will be saved. Ex: "Gut1"
sequencing_depth_histogram = function(DataGG, bins, title, fillcolor, bordercolor, savePlot, saveName){
  DataGGSum<-aggregate(DataGG$GeneCount, by=list(SampleID=DataGG$SampleID), FUN=sum)
  p2 <- ggplot(DataGGSum,aes(x=x)) + 
            geom_histogram(bins = bins, fill=fillcolor, color=bordercolor) + 
            theme(plot.title = element_text(hjust = 0.5)) + 
            labs(title=title, x = "Sequencing depth", y = "Frequency")
  
  print(p2) 
  # Saving plots:
  if(savePlot == TRUE){
    path_save2 <- str_c("../../Result/", saveName,"/", saveName, "_sequencing_depth_histogram.pdf")
    ggsave(filename = path_save2, plot = p2, height = 5, width = 6)
    dev.off()
    print(p2) 
  }
}


## NOT DONE ##

# Performs the entire exploration of a given dataset
# Input: Data = the dataset to be explored (rownames required)
# Output: rData = a logical vector with 1 if a gene should be removed and 0 if it should be kept (regarding low counts) 
explore_dataset = function(Data){
  
  DataGG <- reshaping_dataset(Data)
  
  rData <- compute_low_counts(Data)
  
  return(rData)
  
}



# ==================================== Loading Datasets ====================================#
Gut1<-read.table("../../Data/Raw_data/HumanGutI_COGcountsRaw.txt",header = T,row.names = 1)
Gut2<-read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt",header = T,row.names = 1)
Marine<-read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt",header = T,row.names = 1)

# ===================================== Colour selection ===================================#
# Default colors
colorScale9<-c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
color1<-colorScale9[7]
colorScale150<-colorRampPalette(brewer.pal(9, "YlGnBu"))(150)
colorScale250<-colorRampPalette(brewer.pal(9, "YlGnBu"))(250)
colorScale500<-colorRampPalette(brewer.pal(9, "YlGnBu"))(500)

#==================================== Exploration of data ==================================#

# Reshaping datasets for ggplot
Gut1GG<-reshaping_dataset(Gut1)
Gut2GG<-reshaping_dataset(Gut2)
MarineGG<-reshaping_dataset(Marine)

# Computing low counts
rGut1<-compute_low_counts(Gut1)
rGut2<-compute_low_counts(Gut2)
rMarine<-compute_low_counts(Marine)

# check ordered sequencing depth for Gut1 
Gut1Depth=data.frame(colnames(Gut1), as.numeric(colSums(Gut1)))
gut1Order <- order(Gut1Depth[,2])
Gut1DepthSorted <- Gut1Depth[gut1Order,]  

### For all datasets
# General info
numberOfSamples <- c(ncol(Gut1),ncol(Gut2),ncol(Marine))
numberOfGenes <- c(nrow(Gut1),nrow(Gut2),nrow(Marine))
lowCountGenes <- c(sum(rGut1!=0),sum(rGut2!=0),sum(rMarine!=0))

Dimensions <- cbind(numberOfSamples, numberOfGenes, lowCountGenes)
colnames(Dimensions)<-c("Total number of samples", "Total number of genes", "Genes with low counts")
rownames(Dimensions)<-c("Gut I", "Gut II", "Marine")
print(xtable(Dimensions))

# Summary of Sequensing depths (reads per sample)
SeqSummary<-rbind(summary(colSums(Gut1)),summary(colSums(Gut2)),summary(colSums(Marine)))
rownames(SeqSummary)<-c("Human Gut I", "Human Gut II", "Marine")
print(xtable(SeqSummary))

# Summary of reads per gene
GeneSummary<-rbind(summary(rowSums(Gut1)),summary(rowSums(Gut2)),summary(rowSums(Marine)))
rownames(GeneSummary)<-c("Human Gut I", "Human Gut II", "Marine")
print(xtable(GeneSummary))

# Boxplots for genes in each dataset
gene_boxplot(Gut1GG,"Number of reads in Human Gut I", colorScale150, savePlot=TRUE, saveName="Gut1")
gene_boxplot(Gut2GG,"Number of reads in Human Gut II", colorScale150, savePlot=TRUE, saveName="Gut2")
gene_boxplot(MarineGG,"Number of reads in Marine", colorScale250, savePlot=TRUE, saveName="Marine")

# Histogram of Sequensing depth in each dataset
sequencing_depth_histogram(Gut1GG, bins=30, "Sequencing Depth in Human Gut I", 
  color1, "black", savePlot=TRUE, saveName="Gut1")
sequencing_depth_histogram(Gut2GG, bins=30, "Sequencing Depth in Human Gut II", 
  color1, "black", savePlot=TRUE, saveName="Gut2")
sequencing_depth_histogram(MarineGG, bins=30, "Sequencing Depth in Marine", 
  color1, "black", savePlot=TRUE, saveName="Marine")
