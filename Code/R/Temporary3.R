# Scirpt for creating trade-off tables

library(xtable)

# Read AUC-tables  
AUCq15 <- read.csv("../../Result/Marine/AUC_10q15.csv", header = T)[,c("plotMD","AUC1")]
AUCq30 <- read.csv("../../Result/Marine/AUC_10q30.csv", header = T)[,c("plotMD","AUC1")]

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
#-----------------------------------Trade-off tables for strata --------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

# Read AUC-tables  
AUCAbun3M <- read.csv("../../Result/Marine/AUC_Abundance_3000000_10q15.csv", header = T)[,c("plotMD","AUC1","strata")]
AUCVar3M <- read.csv("../../Result/Marine/AUC_Variability_3000000_10q15.csv", header = T)[,c("plotMD","AUC1","strata")]

# ABUNDANCE
AUCAbun3M <- AUCAbun3M[c(which(AUCAbun3M$plotMD=="m=3, d=1 M"), which(AUCAbun3M$plotMD=="m=6, d=500 k"), 
                    which(AUCAbun3M$plotMD=="m=15, d=200 k"), which(AUCAbun3M$plotMD=="m=30, d=100 k")), ]

AUCVar3M <- AUCVar3M[c(which(AUCVar3M$plotMD=="m=3, d=1 M"), which(AUCVar3M$plotMD=="m=6, d=500 k"), 
                       which(AUCVar3M$plotMD=="m=15, d=200 k"), which(AUCVar3M$plotMD=="m=30, d=100 k")), ]

AUCStrata <- cbind(AUCAbun3M[AUCAbun3M$strata==1,c(1,2)],AUCAbun3M[AUCAbun3M$strata==2,2],
                   AUCAbun3M[AUCAbun3M$strata==3,2], AUCVar3M[AUCVar3M$strata==1,2],
                   AUCVar3M[AUCVar3M$strata==2,2],AUCVar3M[AUCVar3M$strata==3,2])
colnames(AUCStrata) <- c("Experimental Design", "Low Abundance", "Medium Abundance", "High Abundance", 
                            "Low Variability", "Medium Variability", "High Variability")

# print LaTex-code for the table
print(xtable(AUCStrata, digits = 3), include.rownames = F)

