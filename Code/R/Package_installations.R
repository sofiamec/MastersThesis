# Packages required to run Run_entire_analysis.R

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

# For computing AUCs with trapeziodal rule
install.packages("pracma")

# for loading excel-files
install.packages("readxl")

# Chack! Cannot remember
install.packages("plyr")
