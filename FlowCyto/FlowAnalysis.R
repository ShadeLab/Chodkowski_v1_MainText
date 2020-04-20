Flowdata <- read.csv(file="FlowCyto/FlowResults.csv", header=TRUE, sep=",")

#Keep columns of interest
flowData <- Flowdata[,c(1,4:6)]

#Average tech reps within a time point
library(dplyr)
detach("package:plyr", unload=TRUE)
groupedFlowData <- flowData %>% group_by(Species,Biorep,TimePoint)
#Get mean and std
groupedAvgCounts<- summarise(groupedFlowData, mean=mean(Live.CFU.mL))

library(ggplot2)
g=ggplot(groupedAvgCounts, aes(x=as.factor(TimePoint), y=log10(mean)))+
  geom_boxplot(aes(x=as.factor(TimePoint), y=log10(mean)))+
  geom_jitter(aes(x=as.factor(TimePoint), y=log10(mean)))+
  facet_grid(.~Species,)+
  theme_bw(base_size=10)+
  scale_fill_manual(values = "red") + 
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),strip.text.x = element_text(size = 18),strip.text = element_text(face = "italic"))
flowPlot <- g + ylim(c(7,10)) + xlab("Time (h)") + ylab("Log CFU/ml")

ggsave("Figures/SFig4/CellCounts.eps",plot= flowPlot,device="eps",width=15,height=10,dpi=600)


