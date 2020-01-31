# Exploration of data

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(xtable)

# Loading Datasets
gut1<-read.table("MastersThesis/Data/Raw_data/HumanGutI_COGcountsRaw.txt",header = T,row.names = 1)
gut2<-read.table("MastersThesis/Data/Raw_data/HumanGutII_COGcountsRaw.txt",header = T,row.names = 1)
marine<-read.table("MastersThesis/Data/Raw_data/Marine_COGcountsRaw.txt",header = T,row.names = 1)

# Reshaping datasets for ggplot
gut1GG<-read.table("MastersThesis/Data/Raw_data/HumanGutI_COGcountsRaw.txt",header = T)
names(gut1GG)[1] <- "Gene"
gut1GG<-melt(gut1GG, id=c("Gene"))
colnames(gut1GG)<-c("Gene","SampleID","GeneCount")

gut2GG<-read.table("MastersThesis/Data/Raw_data/HumanGutII_COGcountsRaw.txt",header = T)
names(gut2GG)[1] <- "Gene"
gut2GG<-melt(gut2GG, id=c("Gene"))
colnames(gut2GG)<-c("Gene","SampleID","GeneCount")

marineGG<-read.table("MastersThesis/Data/Raw_data/Marine_COGcountsRaw.txt",header = T)
names(marineGG)[1] <- "Gene"
marineGG<-melt(marineGG, id=c("Gene"))
colnames(marineGG)<-c("Gene","SampleID","GeneCount")

# Default colors
colorScale9<-c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
color1<-colorScale9[7]
colorScale150<-colorRampPalette(brewer.pal(9, "YlGnBu"))(150)
colorScale250<-colorRampPalette(brewer.pal(9, "YlGnBu"))(250)
colorScale500<-colorRampPalette(brewer.pal(9, "YlGnBu"))(500)

#===================================================== Exploration of data ==================================================#

# Summary of low counts

## FUNCTION ##
# For a given dataset, this function computes genes which will be removed when filtering.
# It also prints how many genes will be removed respectively kept.
# input: data = the data to be filtered
# output: r = a logical vector with 1 if a gene should be removed, and 0 if it should be kept. 
Compute_Low_Counts=function(data){
  
  a=rowSums(data)<3
  b=vector()
  r=vector()
  for (i in 1:nrow(data)) {
    b[i]<-sum(data[i,]==0)/ncol(data)>=0.75
    r[i]<-a[i]+b[i]
  }
  cat(sprintf("Removed: %s\n", sum(r!=0)))
  cat(sprintf("Kept: %s\n", sum(r==0)))
  return(r)
}

rGut1<-Compute_Low_Counts(gut1)
rGut2<-Compute_Low_Counts(gut2)
rMarine<-Compute_Low_Counts(marine)

# General info
NumberOfSamples<-c(ncol(gut1),ncol(gut2),ncol(marine))
NumberOfGenes<-c(nrow(gut1),nrow(gut2),nrow(marine))
LowCountGenes<-c(sum(rGut1!=0),sum(rGut2!=0),sum(rMarine!=0))
Dimensions<-cbind(NumberOfSamples, NumberOfGenes, LowCountGenes)
colnames(Dimensions)<-c("Total number of samples", "Total number of genes", "Genes with low counts")
rownames(Dimensions)<-c("Gut I", "Gut II", "Marine")
print(xtable(Dimensions))

# Summary of Sequensing depths (reads per sample)
SeqSummary<-rbind(summary(colSums(gut1)),summary(colSums(gut2)),summary(colSums(marine)))
rownames(SeqSummary)<-c("Gut I", "Gut II", "Marine")
print(xtable(SeqSummary))

# Summary of reads per gene
GeneSummary<-rbind(summary(rowSums(gut1)),summary(rowSums(gut2)),summary(rowSums(marine)))
rownames(GeneSummary)<-c("Gut I", "Gut II", "Marine")
print(xtable(GeneSummary))



# Boxplots for genes in Human Gut I
ggplot(gut1GG,aes(x=SampleID,y=GeneCount)) + 
  geom_boxplot(aes(fill=SampleID)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),legend.position = "none") + 
  labs(title="Number of reads in Human Gut I") + 
  scale_fill_manual(values=colorScale150)

# Boxplots for genes in Human Gut II
ggplot(gut2GG,aes(x=SampleID,y=GeneCount)) + 
  geom_boxplot(aes(fill=SampleID)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),legend.position = "none") + 
  labs(title="Number of reads in Human Gut II") + 
  scale_fill_manual(values=colorScale150)

# Boxplots for genes in Marine
ggplot(marineGG,aes(x=SampleID,y=GeneCount)) + 
  geom_boxplot(aes(fill=SampleID)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),legend.position = "none") + 
  labs(title="Number of reads in Marine") + 
  scale_fill_manual(values=colorScale250)


# Histogram of Sequensing depth in Human Gut I
gut1GGSum1<-aggregate(gut1GG$GeneCount, by=list(SampleID=gut1GG$SampleID), FUN=sum)
ggplot(gut1GGSum1,aes(x=x)) + 
  geom_histogram(bins = 60, fill=color1, color="black") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title="Sequencing Depth in Human Gut I", x = "Sequencing depth", y = "Frequency")

# Histogram of Sequensing depth in Human Gut II
gut2GGSum1<-aggregate(gut2GG$GeneCount, by=list(SampleID=gut2GG$SampleID), FUN=sum)
ggplot(gut2GGSum1,aes(x=x)) + 
  geom_histogram(bins = 30, fill=color1, color="black") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title="Sequencing Depth in Human Gut II", x = "Sequencing depth", y = "Frequency")

# Histogram of Sequensing depth in Human Gut I
marineGGSum1<-aggregate(marineGG$GeneCount, by=list(SampleID=marineGG$SampleID), FUN=sum)
ggplot(marineGGSum1,aes(x=x)) +
  geom_histogram(bins = 30, fill=color1, color="black") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title="Sequencing Depth in Marine", x = "Sequencing depth", y = "Frequency")


