
# for plotting
install.packages("ggplot2")

# for reshaping the dataset to works for ggplot
install.packages("reshape2")

# for visualisation in nice colours
install.packages("RColorBrewer")
install.packages("viridis")

# for printing tables compatible with LaTeX
install.packages("xtable")

# enables the use of strings (for plot titles etc.)
install.packages("tidyverse")

# for analysing the datasets and identify DAGs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# for analysing the datasets and identify DAGs
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

# for calculating overdispersion in edgeR
#install.packages("statmod")

# For computing AUCs with trapeziodal rule
install.packages("pracma")
