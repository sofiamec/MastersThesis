#------------------------------------------- Description of Script -----------------------------------------------#

# Scirpt for creating example plots as well as trade-off tables (regular + strata) and cost results for Marine and Gut2 (with DESeq2)

#------------------------------------------------- Librariess ----------------------------------------------------#

library(xtable)
library(ggplot2)
library(viridis)
library(pracma)

#------------------------------------------------- Functions -----------------------------------------------------#

## FUNCTION for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}


## FUNCTION to generate trade-off AUC tables for strata
# input:  reads = total amount of reads (3000000 or 5000000)
#         q = q-value (1.5 or 3)
# output: table containing the trade-off designs and their AUC0.01 values for each stratum

AUC_strata=function(reads,q){
  
  # Read AUC-tables  
  AUCAbun <- read.csv(sprintf("../../Result/%s_DESeq/AUC_Abundance_%d_10q%d.csv", saveName, reads,q*10), header = T)[,c("plotMD","AUC1","strata")]
  AUCVar <- read.csv(sprintf("../../Result/%s_DESeq/AUC_Variability_%d_10q%d.csv",saveName, reads,q*10), header = T)[,c("plotMD","AUC1","strata")]
  
  if(reads==3000000){
    AUCAbun <- AUCAbun[c(which(AUCAbun$plotMD=="m=3, d=1 M"), which(AUCAbun$plotMD=="m=6, d=500 k"), 
                         which(AUCAbun$plotMD=="m=15, d=200 k"), which(AUCAbun$plotMD=="m=30, d=100 k")), ]
    
    AUCVar <- AUCVar[c(which(AUCVar$plotMD=="m=3, d=1 M"), which(AUCVar$plotMD=="m=6, d=500 k"), 
                       which(AUCVar$plotMD=="m=15, d=200 k"), which(AUCVar$plotMD=="m=30, d=100 k")), ]
  } else {
    AUCAbun <- AUCAbun[c(which(AUCAbun$plotMD=="m=5, d=1 M"), which(AUCAbun$plotMD=="m=10, d=500 k"), 
                         which(AUCAbun$plotMD=="m=20, d=250 k"), which(AUCAbun$plotMD=="m=50, d=100 k")), ]
    
    AUCVar <- AUCVar[c(which(AUCVar$plotMD=="m=5, d=1 M"), which(AUCVar$plotMD=="m=10, d=500 k"), 
                       which(AUCVar$plotMD=="m=20, d=250 k"), which(AUCVar$plotMD=="m=50, d=100 k")), ]
  }
  
  AUCTable <- cbind(AUCAbun[AUCAbun$strata==1,c(1,2)],AUCAbun[AUCAbun$strata==2,2],
                    AUCAbun[AUCAbun$strata==3,2], AUCVar[AUCVar$strata==1,2],
                    AUCVar[AUCVar$strata==2,2],AUCVar[AUCVar$strata==3,2])
  colnames(AUCTable) <- c("Experimental Design", "Low Abundance", "Medium Abundance", "High Abundance", 
                          "Low Variability", "Medium Variability", "High Variability")
  
  # print LaTex-code for the table
  cat(sprintf("======== AUC table for %s strata, tradeoff = %d q = %g  =========\n", saveName, reads, q))
  print(xtable(AUCTable, digits = 3), include.rownames = F)
  cat("\n")
  
  return(AUCTable)
}

## FUNCTION for plotting AUC per cost for different designs
AUCPerCost_barplot <- function(Data, tradeoff){
  if (saveName=="Marine"){
    title = bquote("Performance / $"*1000 * " for Marine") # tidigare " $ for Marine"
    yMax = 0.20
  } else if (saveName=="Gut2"){
    title = bquote("Performance / $"*1000 * " for Human Gut II")
    yMax = 0.25
  }
  AUC_cost_plot <- ggplot() + geom_bar(mapping = aes(x=Data$plotMD, y=Data$costAUC, group = Data$q, fill = Data$q), stat = "identity", position=position_dodge()) +
    labs( title = title, fill = "Effect",  x = "Experimental design", y = bquote("Performance / $"*1000)) + # tidigare / "*10^3*" $"
    scale_fill_viridis_d(begin = 0.8, end = 0.6) + ylim(0, yMax)
  
  ggsave(filename = sprintf("../../Result/%s_DESeq/AUCPerCostplot_%s.pdf", saveName, tradeoff), plot = AUC_cost_plot, height = 5, width = 6)
  print(AUC_cost_plot)
}


#-----------------------------------------------------------------------------------------------------------------#
#--------------------------------------- ROC and AUC Example plots -----------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

## Creating three cases for example ROC-plots

#Bad
TPR<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45)/45
FPR<-c(0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45)/45
Class<-rep("Bad (Random)",90)
ROCData1<-data.frame(TPR,FPR,Class)  

#Perfect
TPR<-c(1:30, rep(30,30))/30
FPR<-c(rep(0,30),1:30)/30
Class<-rep("Perfect",60)
ROCData3<-data.frame(TPR,FPR,Class)

#Good
TPR<-c(1,2,3,4,5,5,6,7,8,8,9,10,11,12,12,13,14,15,15,16,17,17,18,18,19,20,21,21,22,23,24,25,25,26,26,27,28,28,28,29,29,29,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,rep(33,30))/33
FPR<-c(0,0,0,0,0,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,6,6,7,7,7,7,7,8,8,9,9,9,10,11,11,12,13,13,14,15,16,17,18,19,19,20,21,22,23,24,25,26,27,28,29,30,31,31,32,33,34,35,36,37,38,38:67)/67
Class<-rep("Good",100)
ROCData2<-data.frame(TPR,FPR,Class)

## Creating Example AUC-plot based on the "Good" ROC-curve
ROCData2_tot <- ROCData2
ROCData2_1<-data.frame(TPR[1:6],FPR[1:6],"1")
colnames(ROCData2_1)<-colnames(ROCData2)
ROCData2_all<-rbind(ROCData2_tot,ROCData2_1)

# Computing AUC values
AUCtot<-trapz(FPR,TPR)
AUC1<-trapz(FPR[1:6],TPR[1:6])/0.01

# Print AUC plot
AUCplot <- ggplot(data=ROCData2_all, aes(x=FPR, y=TPR, group=Class, fill=Class)) +  geom_line(color='black') + 
  geom_ribbon(aes(ymin=0, ymax=TPR), alpha=0.5) + 
  scale_fill_viridis(begin = 0.35, end = 1, discrete=TRUE) +
  labs(x = "False Positive Rate", y = "True Positive Rate", fill = "AUC values", tag = bquote(atop("AUC"[tot]*"   = 0.900", " AUC"[0.01]*" =  0.002"))) +
  theme(legend.text=element_text(color="white"), plot.tag = element_text(size=10), plot.tag.position = c(0.987, 0.513), plot.margin=margin(t = 0.2, r = 1.5, b = 0.2, l = 0.2, unit = "cm")) + guides(color = FALSE) +
  ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))
print(AUCplot)
ggsave(filename = "../../Result/Example_AUCplot.pdf", plot = AUCplot, height = 5, width = 7)

## Creating Example ROC-plot based on the "Good" "Bad" and "Perfect" ROC-curves
ROCData<-rbind(ROCData1,ROCData2,ROCData3)

ROCplot <- ggplot(data=ROCData, aes(x=FPR, y=TPR, group=Class)) +  geom_line(aes(color=Class)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis(begin = 0, end = 0.85, discrete=TRUE) + scale_fill_viridis(begin = 0, end = 0.85, discrete=TRUE) +
  labs(colour=sprintf("Classifier Performance"), x = "False Positive Rate", y = "True Positive Rate") +
  ylim(0, 1) + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) #+ theme(legend.position = c(0.82, 0.18))
print(ROCplot)

ggsave(filename = "../../Result/Example_ROCplot.pdf", plot = ROCplot, height = 5, width = 7)

#-----------------------------------------------------------------------------------------------------------------#
#------------------------------------------- Dataset selection ---------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

Datasets = c("Marine", "Gut2")

for (saveName in Datasets) {

  
  #-----------------------------------------------------------------------------------------------------------------#
  #---------------------------------------------- Trade-off tables -------------------------------------------------#
  #-----------------------------------------------------------------------------------------------------------------#
  
  # Read AUC-tables 
  AUCq15 <- read.csv(sprintf("../../Result/%s_DESeq/AUC_10q15.csv", saveName), header = T)[,c("plotMD","AUC1")]
  AUCq30 <- read.csv(sprintf("../../Result/%s_DESeq/AUC_10q30.csv", saveName), header = T)[,c("plotMD","AUC1")]
  
  # Create table for trade-off for q=1.5
  AUC15 <- cbind(AUCq15[c(which(AUCq15$plotMD=="m=3, d=1 M"), which(AUCq15$plotMD=="m=6, d=500 k"), 
                          which(AUCq15$plotMD=="m=15, d=200 k"), which(AUCq15$plotMD=="m=30, d=100 k")), ],
                 AUCq15[c(which(AUCq15$plotMD=="m=5, d=1 M"), which(AUCq15$plotMD=="m=10, d=500 k"), 
                          which(AUCq15$plotMD=="m=20, d=250 k"), which(AUCq15$plotMD=="m=50, d=100 k")), ])
  colnames(AUC15) <-c("Experimental Desgin (3M)", "AUC0.01", "Experimental Design (5M)","AUC0.01")
  
  # Create table for trade-off for q=3
  AUC3 <- cbind(AUCq30[c(which(AUCq30$plotMD=="m=3, d=1 M"), which(AUCq30$plotMD=="m=6, d=500 k"), 
                         which(AUCq30$plotMD=="m=15, d=200 k"), which(AUCq30$plotMD=="m=30, d=100 k")), ],
                AUCq30[c(which(AUCq30$plotMD=="m=5, d=1 M"), which(AUCq30$plotMD=="m=10, d=500 k"), 
                         which(AUCq30$plotMD=="m=20, d=250 k"), which(AUCq30$plotMD=="m=50, d=100 k")), ])
  colnames(AUC3) <-c("Experimental Desgin (3M)", "AUC0.01", "Experimental Design (5M)","AUC0.01")
  
  # print LaTex-code for the tables
  cat(sprintf("======== AUC tradeoff-table for %s q = 1.5 =========\n", saveName))
  print(xtable(AUC15, digits = 2), include.rownames = F)
  cat("\n")
  cat(sprintf("======== AUC tradeoff-table for %s q = 3 =========\n", saveName))
  print(xtable(AUC3, digits = 2), include.rownames = F)
  cat("\n")
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #----------------------------------- Trade-off tables for strata -------------------------------------------------#
  #-----------------------------------------------------------------------------------------------------------------#
  
  # Generate tables:
  AUC3Mq15=AUC_strata(3000000,1.5)
  AUC3Mq30=AUC_strata(3000000,3)
  AUC5Mq15=AUC_strata(5000000,1.5)
  AUC5Mq30=AUC_strata(5000000,3)
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #--------------------------------------- Cost plots and tables ---------------------------------------------------#
  #-----------------------------------------------------------------------------------------------------------------#
  
  # Structuring the data
  Data15<-read.csv(sprintf("../../Result/%s_DESeq/AUC_10q15.csv", saveName))[,-1]
  Data30<-read.csv(sprintf("../../Result/%s_DESeq/AUC_10q30.csv", saveName))[,-1]
  Data<-rbind(data.frame(Data15, MD=Data15$m*Data15$d, q=1.5), data.frame(Data30, MD=Data15$m*Data15$d, q=3))
  cost=(2*Data$m*100 + Data$d/10000*2*Data$m)
  Data<-data.frame(Data, cost, AUCcost = (cost/Data$AUC1/1000), costAUC = (1000*Data$AUC1/cost))
  Data$plotMD <- factor(Data$plotMD, levels = c("m=3, d=10 k", "m=3, d=100 k", "m=3, d=500 k", "m=3, d=1 M", "m=3, d=5 M", "m=3, d=10 M",
                                                "m=5, d=10 k", "m=5, d=100 k", "m=5, d=500 k", "m=5, d=1 M", "m=5, d=5 M", "m=5, d=10 M",
                                                "m=6, d=500 k", 
                                                "m=10, d=10 k", "m=10, d=100 k", "m=10, d=500 k", "m=10, d=1 M", "m=10, d=5 M", "m=10, d=10 M",  
                                                "m=15, d=200 k", "m=20, d=250 k",     
                                                "m=30, d=10 k", "m=30, d=100 k", "m=30, d=500 k", "m=30, d=1 M", "m=30, d=5 M", "m=30, d=10 M",
                                                "m=50, d=10 k", "m=50, d=100 k", "m=50, d=500 k", "m=50, d=1 M", "m=50, d=5 M", "m=50, d=10 M"))
  Data$q<-as.factor(Data$q)
  Data3M <- Data[Data$MD==3000000,] 
  Data5M <- Data[Data$MD==5000000,] 
  
  # Plotting AUC per Cost barplots
  AUCPerCost_barplot(Data3M, "3M")
  AUCPerCost_barplot(Data5M, "5M")
  
  # Formatting data for tables
  Data$d<-as.factor(Data$d)
  Data$m<-as.factor(Data$m)
  Data15<-Data[Data$q==1.5,]
  Data30<-Data[Data$q==3,]
  
  # Print tables for trade-off 3M
  Data3M$AUC1<-as.character(round2(Data3M$AUC1,2))
  Data3M$cost<-as.character(Data3M$cost)
  Data3M <- Data3M[order(Data3M$plotMD),]
  table3M<-data.frame(Data3M[Data3M$q==1.5,c(9,12,1,14)],Data3M[Data3M$q==3,c(1,14)])
  colnames(table3M)<-c("Experimental design", "Sequencing cost", "AUC$_{0.01}$, q=1.5", "Performance per 1000 US dollar, q=1.5", "AUC$_{0.01}$, q=3", "Performance per 1000 US dollar, q=3")
  
  cat(sprintf("======== Cost table for %s tradeoff 3M =========\n", saveName))
  print(xtable(table3M, digits = 2), include.rownames = F)
  cat("\n")
  
  # Print tables for trade-off 5M
  Data5M$AUC1<-as.character(round2(Data5M$AUC1,2))
  Data5M$cost<-as.character(Data5M$cost)
  Data5M <- Data5M[order(Data5M$plotMD),]
  table5M<-data.frame(Data5M[Data5M$q==1.5,c(9,12,1,14)],Data5M[Data5M$q==3,c(1,14)])
  colnames(table5M)<-colnames(table3M)
  
  cat(sprintf("======== Cost table for %s tradeoff 5M =========\n", saveName))
  print(xtable(table5M, digits = 2), include.rownames = F)
  cat("\n")

}
