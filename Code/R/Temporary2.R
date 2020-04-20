library(ggplot2)

Data15<-read.csv("../../Result/Marine_DESeq_Median_100Rep/AUC_10q15.csv")[,-1]

Data30<-read.csv("../../Result/Marine_DESeq_Median_100Rep/AUC_10q30.csv")[,-1]

Data15<-data.frame(Data15, MD=Data15$m*Data15$d,q=1.5)
Data30<-data.frame(Data30, MD=Data15$m*Data15$d,q=3)
Data<-rbind(Data15,Data30)

cost=(2*Data$m*100+Data$d/10000*2*Data$m)
Data<-data.frame(Data, cost, AUCcost = (cost/Data$AUC1))

Data$plotMD <- factor(Data$plotMD, levels = c("m=3 d=10k", "m=3 d=100k", "m=3 d=500k", "m=3 d=1M", "m=3 d=5M", "m=3 d=10M",
                                              "m=5 d=10k", "m=5 d=100k", "m=5 d=500k", "m=5 d=1M", "m=5 d=5M", "m=5 d=10M",
                                              "m=6 d=500k", 
                                              "m=10 d=10k", "m=10 d=100k", "m=10 d=500k", "m=10 d=1M", "m=10 d=5M", "m=10 d=10M",  
                                              "m=15 d=200k", "m=20 d=250k",     
                                              "m=30 d=10k", "m=30 d=100k", "m=30 d=500k", "m=30 d=1M", "m=30 d=5M", "m=30 d=10M",
                                              "m=50 d=10k", "m=50 d=100k", "m=50 d=500k", "m=50 d=1M", "m=50 d=5M", "m=50 d=10M"))

Data3M <- Data[Data$MD==3000000,] 
Data5M <- Data[Data$MD==5000000,] 

ggplot() + geom_bar(mapping = aes(x=Data3M$plotMD, y=Data3M$AUC1, group = Data3M$q, fill = Data3M$q), stat = "identity", position=position_dodge()) + 
  geom_line(mapping = aes(x=Data3M$plotMD, y=Data3M$cost/6000, group = 1), size = 1, color = "blue") +
  scale_y_continuous(name = "AUC1", 
                     sec.axis = sec_axis(~.*6000, name = "Sequencing cost"))+#, 
                                         #labels = function(b) { paste0(round(b * 100, 0), "%")})) + 
  theme(
    axis.title.y = element_text(color = "grey"),
    axis.title.y.right = element_text(color = "blue"))

ggplot() + geom_line(mapping = aes(x=Data3M$plotMD, y=Data3M$cost/1, group = 1), size = 1, color = "blue")

ggplot(Data3M[Data3M$q==1.5,]) + geom_line(aes(plotMD, cost, group=1))
