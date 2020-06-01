# Packages required to run Exploration_of_data.R and Master.R

# For plotting
install.packages("ggplot2")

# For reshaping the dataset to works for ggplot
install.packages("reshape2")

# For visualisation in nice colours
install.packages("RColorBrewer")
install.packages("viridis")

# For printing tables compatible with LaTeX
install.packages("xtable")

# Enables the use of strings (for plot titles etc.)
install.packages("tidyverse")

# For analysing the datasets and identify DAGs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# For computing AUCs with trapeziodal rule
install.packages("pracma")

# For loading excel-files
install.packages("readxl")

# For handling ROC and AUC matrixes
install.packages("plyr")
