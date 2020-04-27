library(ggplot2)
library(viridis)
library(xtable)

## FUNCTION for always rounding 0.5, 0.05 etc upwards
round2 <- function(x, n) {
  z = trunc(abs(x)*10^n +0.5)/10^n *sign(x)
  return(z)
}

## FUNCTION for plotting cost vs performance in a combined bar- and lineplot
cost_barplot <- function(Data, tradeoff){
  costplot <-ggplot() + geom_bar(mapping = aes(x=Data$plotMD, y=Data$AUC1, group = Data$q, fill = Data$q), stat = "identity", position=position_dodge()) + 
    geom_line(mapping = aes(x=Data$plotMD, y=Data$cost/max(Data$cost), group = 1), size = 1, color = "#9E0142") +
    scale_y_continuous(name = "Performance in AUC1\n", sec.axis = sec_axis(~.*max(Data$cost), name = "Sequencing cost in US dollars\n")) + 
    theme( axis.title.y = element_text(color = "#00BC77"),
           axis.title.y.right = element_text(color = "#9E0142"),
           axis.text.y = element_text(color="#00BC77"),
           axis.text.y.right = element_text(color="#9E0142"),) +
    labs( title = bquote("Cost vs. performance in"~AUC[0.01]), fill = "q-value",  x = "Experimental design") + 
    scale_fill_viridis_d(begin = 0.8, end = 0.7)
  
  ggsave(filename = sprintf("../../Result/Costplot_Marine%s.pdf", tradeoff), plot = costplot, height = 5, width = 6)
  print(costplot)
}

## FUNCTION for plotting cost per AUC for different designs
CostPerAUC_barplot <- function(Data, tradeoff){
  cost_AUC_plot <- ggplot() + geom_line(mapping = aes(x=Data$plotMD, y=Data$AUCcost, group = Data$q, colour = Data$q), size = 1) +
    labs( title = bquote("Cost per AUC based on"~AUC[0.01]), fill = "q-value",  x = "Experimental design") + 
    scale_fill_viridis_d(begin = 0.8, end = 0.7) 
  
  ggsave(filename = sprintf("../../Result/CosPerAUCtplot_Marine%s.pdf", tradeoff), plot = cost_AUC_plot, height = 5, width = 6)
  print(cost_AUC_plot)
}

## FUNCTION for plotting AUC per cost for different designs
AUCPerCost_barplot <- function(Data, tradeoff){
  AUC_cost_plot <- ggplot() + geom_line(mapping = aes(x=Data$plotMD, y=Data$costAUC, group = Data$q, colour = Data$q), size = 1) +
    labs( title = bquote("AUC per cost based on"~AUC[0.01]), fill = "q-value",  x = "Experimental design") + 
    scale_fill_viridis_d(begin = 0.8, end = 0.7) 
  
  ggsave(filename = sprintf("../../Result/AUCPerCostplot_Marine%s.pdf", tradeoff), plot = AUC_cost_plot, height = 5, width = 6)
  print(AUC_cost_plot)
}


#############

Data15<-read.csv("../../Result/Temporary/Marine_DESeq/AUC_10q15.csv")[,-1]

Data30<-read.csv("../../Result/Temporary/Marine_DESeq/AUC_10q30.csv")[,-1]

Data15<-data.frame(Data15, MD=Data15$m*Data15$d,q=1.5)
Data30<-data.frame(Data30, MD=Data15$m*Data15$d,q=3)
Data<-rbind(Data15,Data30)

cost=(2*Data$m*100+Data$d/10000*2*Data$m)
Data<-data.frame(Data, cost, AUCcost = (cost/Data$AUC1), costAUC = (Data$AUC1/cost))

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

cost_barplot(Data3M, "3M")
cost_barplot(Data5M, "5M")

CostPerAUC_barplot(Data3M, "3M")
CostPerAUC_barplot(Data5M, "5M")


AUCPerCost_barplot(Data3M, "3M")
AUCPerCost_barplot(Data5M, "5M")

ggplot() + geom_line(mapping = aes(x=Data$plotMD, y=Data$AUCcost, group = Data$q, colour = Data$q), size = 1) +
  labs( title = bquote("Cost vs. performance in"~AUC[0.01]), fill = "q-value",  x = "Experimental design") + 
  scale_fill_viridis_d(begin = 0.8, end = 0.7)

# Tables
Data3M$AUC1<-as.character(round2(Data3M$AUC1,2))
Data3M <- Data3M[order(Data3M$AUCcost),]
Data3M <- Data3M[order(Data3M$q),]
table15 <- Data3M[,-c(2,3,4,5,6,7,8,10)]
print(xtable(table15, digits = 0), include.rownames = F)

Data5M$AUC1<-as.character(round2(Data5M$AUC1,2))
Data5M <- Data5M[order(Data5M$AUCcost),]
Data5M <- Data5M[order(Data5M$q),]
table30 <- Data5M[,-c(2,3,4,5,6,7,8,10)]
print(xtable(table30, digits = 0), include.rownames = F)
