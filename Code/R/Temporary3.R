# Scirpt for creating trade-off tables

library(xtable)

# Read AUC-tables  
AUCq15 <- read.csv("../../Result/Marine_DESeq/AUC_10q15.csv", header = T)[,c("plotMD","AUC1")]
AUCq30 <- read.csv("../../Result/Marine_DESeq/AUC_10q30.csv", header = T)[,c("plotMD","AUC1")]

# Create table for trade-off for 3M reads in total
AUC3M <- cbind(AUCq15[c(which(AUCq15$plotMD=="m=3, d=1 M"), which(AUCq15$plotMD=="m=6, d=500 k"), 
                   which(AUCq15$plotMD=="m=15, d=200 k"), which(AUCq15$plotMD=="m=30, d=100 k")), ],
               AUCq30[c(which(AUCq30$plotMD=="m=3, d=1 M"), which(AUCq30$plotMD=="m=6, d=500 k"), 
                        which(AUCq30$plotMD=="m=15, d=200 k"), which(AUCq30$plotMD=="m=30, d=100 k")), "AUC1"])
colnames(AUC3M) <-c("Experimental Desgin", "AUC0.01 for q=1.5", "AUC0.01 for q=3")

# Create table for trade-off for 5M reads in total
AUC5M <- cbind(AUCq15[c(which(AUCq15$plotMD=="m=5, d=1 M"), which(AUCq15$plotMD=="m=10, d=500 k"), 
                        which(AUCq15$plotMD=="m=20, d=250 k"), which(AUCq15$plotMD=="m=50, d=100 k")), ],
               AUCq30[c(which(AUCq30$plotMD=="m=5, d=1 M"), which(AUCq30$plotMD=="m=10, d=500 k"), 
                        which(AUCq30$plotMD=="m=20, d=250 k"), which(AUCq30$plotMD=="m=50, d=100 k")), "AUC1"])
colnames(AUC5M) <-c("Experimental Desgin", "AUC0.01 for q=1.5", "AUC0.01 for q=3")

# print LaTex-code for the tables
print(xtable(AUC3M, digits = 3), include.rownames = F)
print(xtable(AUC5M, digits = 3), include.rownames = F)

#-----------------------------------------------------------------------------------------------------------------#
#----------------------------------- Trade-off tables for strata -------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Choose which table to generate:
depth=5000000                   # "3000000" or "5000000"
q=30                            # "15"      or "3"

# Read AUC-tables  
AUCAbun <- read.csv(sprintf("../../Result/Marine_DESeq/AUC_Abundance_%d_10q%d.csv",depth,q), header = T)[,c("plotMD","AUC1","strata")]
AUCVar <- read.csv(sprintf("../../Result/Marine_DESeq/AUC_Variability_%d_10q%d.csv",depth,q), header = T)[,c("plotMD","AUC1","strata")]

if(depth==3000000){
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


AUCStrata <- cbind(AUCAbun[AUCAbun$strata==1,c(1,2)],AUCAbun[AUCAbun$strata==2,2],
                   AUCAbun[AUCAbun$strata==3,2], AUCVar[AUCVar$strata==1,2],
                   AUCVar[AUCVar$strata==2,2],AUCVar[AUCVar$strata==3,2])
colnames(AUCStrata) <- c("Experimental Design", "Low Abundance", "Medium Abundance", "High Abundance", 
                            "Low Variability", "Medium Variability", "High Variability")

# print LaTex-code for the table
print(xtable(AUCStrata, digits = 3), include.rownames = F)

