# Scirpt for creating trade-off tables

library(xtable)

# Read AUC-tables q=1.5 
AUCq15 <- read.csv("../../Result/Marine_DESeq/AUC_10q15.csv", header = T)[,c("plotMD","AUC1")]
AUCq30 <- read.csv("../../Result/Marine_DESeq/AUC_10q30.csv", header = T)[,c("plotMD","AUC1")]

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
print(xtable(AUC15, digits = 3), include.rownames = F)
print(xtable(AUC3, digits = 3), include.rownames = F)

#-----------------------------------------------------------------------------------------------------------------#
#----------------------------------- Trade-off tables for strata -------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

## FUNCTION to generate trade-off AUC tables for strata
# input: total amount of reads (3000000 or 5000000), q-value (1.5 or 3)
# output: table containing the trade-off designs and their AUC0.01 values for each stratum

AUC_strata=function(reads,q){
  
  # Read AUC-tables  
  AUCAbun <- read.csv(sprintf("../../Result/Marine_DESeq/AUC_Abundance_%d_10q%d.csv",reads,q), header = T)[,c("plotMD","AUC1","strata")]
  AUCVar <- read.csv(sprintf("../../Result/Marine_DESeq/AUC_Variability_%d_10q%d.csv",reads,q), header = T)[,c("plotMD","AUC1","strata")]
  
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
  print(xtable(AUCTable, digits = 3), include.rownames = F)
  
  return(AUCTable)
}

# Generate tables:
AUC3Mq15=AUC_strata(3000000,15)
AUC5Mq15=AUC_strata(5000000,15)
AUC3Mq30=AUC_strata(3000000,30)
AUC5Mq30=AUC_strata(5000000,30)


