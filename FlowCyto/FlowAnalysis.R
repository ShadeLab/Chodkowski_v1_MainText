Flowdata <- read.csv(file="FlowCyto/FlowResults.csv", header=TRUE, sep=",")

#Keep columns of interest
flowData <- Flowdata[,c(1,5:8)]

#Average tech reps within a time point
library(dplyr)
detach("package:plyr", unload=TRUE)
groupedFlowData <- flowData %>% group_by(Species,Biorep,TimePoint)
#Get mean and std
groupedAvgCounts.live <- summarise(groupedFlowData, mean=mean(Live.CFU.mL))
groupedAvgCounts.dead <- summarise(groupedFlowData, mean=mean(Dead.Population))

groupedAvgCounts.live$viability <- "Live"
groupedAvgCounts.dead$viability <- "Dead"

groupedAvgCounts <- rbind(groupedAvgCounts.live,groupedAvgCounts.dead)

library(ggplot2)
g=ggplot(groupedAvgCounts, aes(x=as.factor(TimePoint), y=log10(mean),fill=viability))+
  geom_boxplot(aes(x=as.factor(TimePoint), y=log10(mean)))+
  geom_jitter(aes(x=as.factor(TimePoint), y=log10(mean)))+
  facet_grid(.~Species,)+
  theme_bw(base_size=10)+
  scale_fill_manual(values = c("blue","green")) +
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),strip.text.x = element_text(size = 18),strip.text = element_text(face = "italic"))
flowPlot <- g + ylim(c(6,10)) + xlab("Time (h)") + ylab("Log cells/mL") +
  theme(legend.title = element_text(size = 16),legend.text = element_text(size = 12),legend.position = "bottom") +
  labs(fill = "Viability assessment") + guides(fill = guide_legend(reverse = TRUE))

ggsave("Figures/SFig4/CellCounts.eps",plot= flowPlot,device="eps",width=15,height=10,dpi=600)
