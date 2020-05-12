
## Exploration of data


#======================================== Libraries ======================================#
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(xtable)
library(tidyverse)
#library(plyr)
library(viridis)
library(DESeq2)
#library(pracma)
library("readxl")
#======================================== Functions ======================================#

## Remove low counts 
# For a given dataset, this function removes genes with low counts (>75 % or an average count <3).
# input: Data = the data to remove genes from
#        filterAll = if T, both kriteria are checked (75%, <3). If F, only kriteria <3 is checked
# output: 
remove_low_counts=function(Data, filterAll){
  a=rowSums(Data)<3
  
  if (filterAll==T){
    b=vector()
    r=vector()
    for (i in 1:nrow(Data)) {
      b[i]<-sum(Data[i,]==0)/ncol(Data)>0.75
      r[i]<-a[i]+b[i]
    }
    FilteredData=Data[r==0,]
    
  } else {
    FilteredData=Data[a==F,]
  }
  return(FilteredData)
}

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

# gene_histogram
# Displays and saves a figure of how the sequencing depth of samples are distributed in a dataset
# Input: DataGG = the reshaped dataset to be plotted
#        bins = number of partitions for the histogram
#        title = a string of what should be displayed as the figure title. Ex: "Sequencing Depth in Human Gut I"
#        fillcolor = colour code for filling the histogram bars. Ex: colour1
#        bordercolor = colour code for filling the bar borders. Ex: "black"
#        savePlot = TRUE or FALSE, statement indicates if the plot should be saved or not
#        saveName = character with name of the dataset as it will be saved. Ex: "Gut1"
gene_histogram = function(DataGG, bins, title, fillcolor, bordercolor, savePlot, saveName){
  DataGGSum<-aggregate(DataGG$GeneCount, by=list(Gene=DataGG$Gene), FUN=sum)
  p2 <- ggplot(DataGGSum,aes(x=x)) + 
    geom_histogram(bins = bins, fill=fillcolor, color=bordercolor) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title=title, x = "Total Gene Count/Abundance", y = "Frequency")
  
  print(p2) 
  # Saving plots:
  if(savePlot == TRUE){
    path_save2 <- str_c("../../Result/", saveName,"/", saveName, "_gene_count_histogram.pdf")
    ggsave(filename = path_save2, plot = p2, height = 5, width = 6)
    dev.off()
    print(p2) 
  }
}

## FUNCTION for creating strata
# DESeq2-analysis for original data 
# This function uses DESeq2 to estimate the mean count and dispersion for genes in a dataset 
# (after normalising them based on sequencing depth). This is then used for creating strata.
#
# Input:  Data = the data to analyse (Gut2 or Marine)
#         numberOfStrata = the number of strata/levels used when dividing the data
# Output: Result = a dataframe containing the mean values and the estimated dispersion for each gene 
#                  as well as the corresponding abundance- and variability-strata for each gene
DESeq2_for_strata=function(Data,numberOfStrata){
  
  DesignMatrix <- data.frame(group=factor(rep(1,(ncol(Data)))))                        # all samples belong to the same group 
  CountsDataset<-DESeqDataSetFromMatrix(countData=Data,DesignMatrix, design=~1) # combine design matrix and data into a dataset
  ResultDESeq<-suppressMessages(DESeq(CountsDataset))                           # Perform analysis (suppress messages from it) 
  Res=results(ResultDESeq,independentFiltering=FALSE,cooksCutoff=FALSE)         # extract results
  
  Result <- data.frame(rownames(Data), Res$baseMean, dispersions(ResultDESeq))  # dataframe with genes, their mean values (after normalization) and dispersion estimates
  Result <- data.frame(Result[,-1], row.names=Result[,1])                       # put first column (genes) as rowname
  
  # create a vector for assigning smallest to highest strata
  n=numberOfStrata
  repVector<-as.integer(cumsum(c(rep(nrow(Data)/n,n))))
  repVector<-repVector-c(0,repVector)[-(length(repVector)+1)]
  strataVector<-rep(1:n,repVector)
  
  # sort by abundance and assign abundance-strata
  Result<-Result[order(Result[,1]),] 
  Result<-cbind(Result,strataVector)
  
  # sort by variability and assign variability-strata
  Result<-Result[order(Result[,2]),]
  Result<-cbind(Result,strataVector)
  
  # rename columns, factorise strata and sort genes by original order
  colnames(Result)<-c("BaseMean", "Estimated dispersion", "AbundanceStrata", "VariabilityStrata")
  Result$AbundanceStrata<-as.factor(Result$AbundanceStrata)
  Result$VariabilityStrata<-as.factor(Result$VariabilityStrata)
  Result<-Result[rownames(Data),]
  
  rm(n, numberOfStrata, repVector,strataVector)
  return(Result)
}

# strata_summary
# This function calculates summary statistics for the different stratas of a stratified dataset
# Input: Data=The stratified dataset to analyse, where the columns NEEDS to be in the following order:
#             BaseMean, Estimated Dispersion, AbundanceStrata, VariabilityStrata
strata_summary=function(Data, numberOfStrata){
  if (numberOfStrata==3){
    StrataSummary <- format(rbind(summary(Data[Data[,3]==1,1]),summary(Data[Data[,3]==2,1]),
                                  summary(Data[Data[,3]==3,1]), summary(Data[Data[,4]==1,2]),
                                  summary(Data[Data[,4]==2,2]), summary(Data[Data[,4]==3,2])),
                            scientific = F, digits = 2)
    rownames(StrataSummary) <- c("Abundance Strata 1", "Abundance Strata 2", "Abundance Strata 3", 
                                 "Disperstion strata 1", "Disperstion strata 2", "Disperstion strata 3")
  } else if (numberOfStrata==5){
    StrataSummary <- format(rbind(summary(Data[Data[,3]==1,1]),summary(Data[Data[,3]==2,1]),
                                  summary(Data[Data[,3]==3,1]), summary(Data[Data[,3]==4,1]),
                                  summary(Data[Data[,3]==5,1]), summary(Data[Data[,4]==1,2]),
                                  summary(Data[Data[,4]==2,2]), summary(Data[Data[,4]==3,2]),
                                  summary(Data[Data[,4]==4,2]), summary(Data[Data[,4]==5,2])),
                            scientific = F, digits = 2)
    rownames(StrataSummary) <- c("Abundance Strata 1", "Abundance Strata 2", "Abundance Strata 3", "Abundance Strata 4", "Abundance Strata 5", 
                                 "Disperstion strata 1", "Disperstion strata 2", "Disperstion strata 3", "Disperstion strata 4", "Disperstion strata 5")
  }

  return(StrataSummary)
}


# ==================================== Loading Datasets ====================================#
Gut2Original<-read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt",header = T,row.names = 1)
MarineOriginal<-read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt",header = T,row.names = 1)
Gut2OriginalPreGG<-read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt",header = T)
MarineOriginalPreGG<-read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt",header = T)

{ ResistanceOriginal=t(read_excel("../../Data/Raw_data/GENE_QUANTIFICATIONS.raw.xlsx")[,-c(2,4)])
  colnames(ResistanceOriginal) <- ResistanceOriginal[1,]
  # Extracting Human samples:
  ResistanceOriginal<-ResistanceOriginal[,ResistanceOriginal[2,]=="Airways" | ResistanceOriginal[2,]=="Gastrointestinal" | ResistanceOriginal[2,]=="Oral" | ResistanceOriginal[2,]=="Skin" | ResistanceOriginal[2,]=="Urogenital"]
  ResistanceOriginal<-ResistanceOriginal[-2,]
  ResistanceOriginal=data.frame(row.names = row.names(ResistanceOriginal)[-1], apply(ResistanceOriginal[-1,],2,as.integer))
  ResistanceOriginalPreGG=data.frame(row.names(ResistanceOriginal),apply(ResistanceOriginal,2,as.integer))
}

ResistanceIntermediate=ResistanceOriginal[,colSums(ResistanceOriginal)>=10000000] # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
Resistance <- remove_low_counts(ResistanceIntermediate, filterAll=F) # This is the dataset used in analysis where samples and genes with low counts are removed

Gut2Intermediate = Gut2Original[,colSums(Gut2Original)>=5000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
Gut2 <- remove_low_counts(Gut2Intermediate, filterAll=T) # This is the dataset used in analysis where samples and genes with low counts are removed

MarineIntermediate = MarineOriginal[,colSums(MarineOriginal)>=10000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
Marine <- remove_low_counts(MarineIntermediate, filterAll=T) # This is the dataset used in analysis where samples and genes with low counts are removed

## MAYBE REPEAT THE EXPLORATION BELOW FOR Gut2 and Marine??
rm(Gut2Intermediate, MarineIntermediate, ResistanceIntermediate)

# ===================================== Colour selection ===================================#
# Default colors
colorScale9<-c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
color1<-colorScale9[7]
colorScale150<-colorRampPalette(brewer.pal(9, "YlGnBu"))(150)
colorScale250<-colorRampPalette(brewer.pal(9, "YlGnBu"))(250)
colorScale500<-colorRampPalette(brewer.pal(9, "YlGnBu"))(500)

#=============================================================================================================================#
#====================================== Exploration of Data ==================================================================#

#====================================== General Exploration ==================================#

# Reshaping datasets for ggplot
Gut2GGOriginal<-reshaping_dataset(Gut2OriginalPreGG)
MarineGGOriginal<-reshaping_dataset(MarineOriginalPreGG)
ResistanceGGOriginal<-reshaping_dataset(ResistanceOriginalPreGG[-1,]) # The non-resistance gene is removed to better compare the other genes in the dataset

# Computing low counts
rGut2<-compute_low_counts(Gut2Original)
rMarine<-compute_low_counts(MarineOriginal)
rResistance <- compute_low_counts(ResistanceOriginal)

### For all datasets
# General info
numberOfSamples <- c(ncol(Gut2Original),ncol(MarineOriginal),ncol(ResistanceOriginal))
numberOfGenes <- c(nrow(Gut2Original),nrow(MarineOriginal), nrow(ResistanceOriginal))
lowCountGenes <- c(sum(rGut2!=0),sum(rMarine!=0), sum(rResistance!=0))

Dimensions <- cbind(numberOfSamples, numberOfGenes, lowCountGenes)
colnames(Dimensions)<-c("Total number of samples", "Total number of genes", "Genes with low counts")
rownames(Dimensions)<-c("Gut II", "Marine", "Resistance")
print(xtable(Dimensions))

# Summary of Sequensing depths (reads per sample)
SeqSummary<-rbind(summary(colSums(Gut2Original)),summary(colSums(MarineOriginal)), summary(colSums(ResistanceOriginal)))
rownames(SeqSummary)<-c("Human Gut II", "Marine", "Resistance")
print(xtable(SeqSummary))

SeqSummary2<-rbind(summary(colSums(Gut2Original)),summary(colSums(MarineOriginal)),summary(colSums(ResistanceOriginal)),summary(colSums(Gut2)),summary(colSums(Marine)),summary(colSums(Resistance)))
rownames(SeqSummary2)<-c("Human Gut II", "Marine", "Antibiotic Resistance Genes", "Filtered Human Gut II", "Filtered Marine", "Filtered Antibiotic Resistance Genes")

# Summary of reads per gene
GeneSummary<-rbind(summary(rowSums(Gut2Original)), summary(rowSums(MarineOriginal)))
rownames(GeneSummary)<-c("Human Gut II", "Marine")
print(xtable(GeneSummary))

GeneSummary2<-rbind(summary(rowSums(Gut2Original)),summary(rowSums(MarineOriginal)),summary(rowSums(Gut2)),summary(rowSums(Marine)))
rownames(GeneSummary2)<-c("Human Gut II", "Marine","Filtered Human Gut II", "Filtered Marine")

# Boxplots for genes in each dataset
gene_boxplot(Gut2GGOriginal,"Number of reads in Human Gut II", colorScale150, savePlot=F, saveName="Gut2")
gene_boxplot(MarineGGOriginal,"Number of reads in Marine", colorScale250, savePlot=F, saveName="Marine")

# Histogram of Sequensing depth in each dataset
sequencing_depth_histogram(Gut2GGOriginal, bins=30, "Sequencing Depth in Human Gut II", 
                           color1, "black", savePlot=F, saveName="Gut2")
sequencing_depth_histogram(MarineGGOriginal, bins=30, "Sequencing Depth in Marine", 
                           color1, "black", savePlot=F, saveName="Marine")
sequencing_depth_histogram(ResistanceGGOriginal, bins=30, "Sequencing Depth in Antibiotic Resistance Genes", 
                           color1, "black", savePlot=F, saveName="Resistance")


# Histogram of Sequensing depth in each dataset
gene_histogram(Gut2GGOriginal, bins=30, "Gene Abundance in Human Gut II", 
                           color1, "black", savePlot=F, saveName="Gut2")
gene_histogram(MarineGGOriginal, bins=30, "Gene Abundance in Marine", 
                           color1, "black", savePlot=F, saveName="Marine")
gene_histogram(ResistanceGGOriginal, bins=30, "Gene Abundance in Antibiotic Resistance Genes", 
                           color1, "black", savePlot=F, saveName="Resistance")


#==================================== Stratified Datasets ==================================#
Gut2Strata<- DESeq2_for_strata(Gut2,3)
MarineStrata<-DESeq2_for_strata(Marine,3)
ResistanceStrata<-DESeq2_for_strata(Resistance,3)

# Table containing a summary for each strata
Gut2StrataSummary=strata_summary(Gut2Strata,3)
MarineStrataSummary=strata_summary(MarineStrata,3)
ResistanceStrataSummary=strata_summary(ResistanceStrata,3)

#=============================================================================================================================#
#====================================== Example tables for Skype ==================================================================#

AUCAb15<-read.csv("../../Result/Gut2/AUC_Abundance_3000000_10q15.csv")[,-1]
AUCV15<-read.csv("../../Result/Gut2/AUC_Variability_3000000_10q15.csv")[,-1]

AUCAb30<-read.csv("../../Result/Gut2/AUC_Abundance_3000000_10q30.csv")[,-1]
AUCV30<-read.csv("../../Result/Gut2/AUC_Variability_3000000_10q30.csv")[,-1]

print(xtable(AUCAb15[AUCAb15$strata==1,c(2,5)],caption = "High abundance"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCAb15$strata==2,c(2,5)],caption = "Medium abundance"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCAb15$strata==3,c(2,5)],caption = "Low abundance"), include.rownames=FALSE)

print(xtable(AUCAb15[AUCV15$strata==1,c(2,5)],caption = "Low variability"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCV15$strata==2,c(2,5)],caption = "Medium variability"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCV15$strata==3,c(2,5)],caption = "High variability"), include.rownames=FALSE)

print(xtable(AUCAb15[AUCAb30$strata==1,c(2,5)],caption = "High abundance"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCAb30$strata==2,c(2,5)],caption = "Medium abundance"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCAb30$strata==3,c(2,5)],caption = "Low abundance"), include.rownames=FALSE)

print(xtable(AUCAb15[AUCV30$strata==1,c(2,5)],caption = "Low variability"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCV30$strata==2,c(2,5)],caption = "Medium variability"), include.rownames=FALSE)
print(xtable(AUCAb15[AUCV30$strata==3,c(2,5)],caption = "High variability"), include.rownames=FALSE)


Genes15<-read.csv("../../Result/Marine/GenesFDR_10q15.csv")[,-1]
Genes30<-read.csv("../../Result/Marine/GenesFDR_10q30.csv")[,-1]

print(xtable(Genes15),include.rownames = F)
print(xtable(Genes30),include.rownames = F)

Genes15<-read.csv("../../Result/Marine_oGLM/GenesFDR_10q15.csv")[,-1]
Genes30<-read.csv("../../Result/Marine_oGLM/GenesFDR_10q30.csv")[,-1]

print(xtable(Genes15),include.rownames = F)
print(xtable(Genes30),include.rownames = F)
