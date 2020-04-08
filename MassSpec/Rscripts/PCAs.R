#Read in data

###Polar Pos###
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Species"=mDataSplit$X2, "Time"=mDataSplitTime$X1)

#Prepare metadata table for varpart
metaDvP <- metaD
#Add row names
row.names(metaDvP) <- mData

metaD <- data.frame("Group"=paste(mDataSplit$X2,mDataSplitTime$X1))
#Remove column names, remove row 1, and transpose dataframe
#colnames(dataNorm) <- NULL
dataNorm <- dataNorm[-1,]
#Convert characters to numeric
dataNorm <- as.data.frame(lapply(dataNorm, as.numeric))

#transpose
dataNormt <- t(dataNorm)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="euclidean")
#groups <- factor(c(rep(1,16), rep(2,8))
groups <- c(metaD$Group)
#mod <- betadisper(dist.Metab ~ Condition + Time, metaD)
mod <- betadisper(dist.Metab, groups)
#Calculates spatial median, not center of mass
#plot(mod)
#pcoa.Metab <- cmdscale(dist.Metab, eig=TRUE)
#Variance explained should match individual plots from Metaboanalyst
PC1var <- mod$eig[1]/sum(mod$eig) #Checks out
PC2var <- mod$eig[2]/sum(mod$eig) #Checks out

#Extract scores

modScores <- scores(mod)
#Extract centroid
centroids <- as.data.frame(modScores$centroids)
#Extract species and time info
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
#Remove redundancies
condI <- condI[!duplicated(condI),]
#Add species
centroids$Label <- condI$Condition
#Add time
centroids$Time <- condI$Time
#Flip y-axis
centroids$PCoA2 <- centroids$PCoA2*-1


#Extract axes vectors from individual samples
sites <- as.data.frame(modScores$sites)
#add groups
sites$groups <- mod$group


#Calculated PCoA axes SDs for each group
ag <- aggregate(.~groups,sites, function(x) sd=sd(x))
#Std axis 1
sd_axis1 <- ag$PCoA1
#Std axis 2
sd_axis2 <- ag$PCoA2

library(ggplot2)
#library(repr)
#options(repr.plot.width = 5, repr.plot.height = 3)
#Change time from factor to numeric
centroids$Time <- rep(c(12.5,25,30,35,40,45),3)

PCA_PolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% var. explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% var. explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 20),legend.title=element_text(size=20))

PCA_PP_Bt <- PCA_PolarPos + ylim(-55,25) + xlim(-200,-75) + theme(axis.title=element_blank())
PCA_PP_Cv <- PCA_PolarPos + ylim(47.5,60) + xlim(42.5,57) + theme(axis.title=element_blank())
PCA_PP_Ps <- PCA_PolarPos + ylim(-80,-30) + xlim(50,90) + theme(axis.title=element_blank())

library(patchwork)
#PCA_PolarPos_plots <- wrap_plots(PCA_PolarPos,PCA_PP_Cv,PCA_PP_Bt,PCA_PP_Ps) +
# theme(legend.position="right",text = element_text(size=20)) + guides(fill=guide_legend(title = "Strain",order=1),size=guide_legend(title = "Time",order=2),color=guide_legend(order=3))

PCA_PolarPos_plots <- PCA_PolarPos + PCA_PP_Cv +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_PP_Bt + PCA_PP_Ps + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')

ggsave("Figures/Fig1/PCA_PolarPos.eps",plot=PCA_PolarPos_plots,device="eps",width=15,height=10,dpi=600)


#Run variation partitioning.
PPvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.63722 #Time: 0.00587 #Species + Time:0.69384





#####Polar Neg#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Species"=mDataSplit$X2, "Time"=mDataSplitTime$X1)

#Prepare metadata table for varpart
metaDvP <- metaD
#Add row names
row.names(metaDvP) <- mData

metaD <- data.frame("Group"=paste(mDataSplit$X2,mDataSplitTime$X1))
#Remove column names, remove row 1, and transpose dataframe
#colnames(dataNorm) <- NULL
dataNorm <- dataNorm[-1,]
#Convert characters to numeric
dataNorm <- as.data.frame(lapply(dataNorm, as.numeric))

#transpose
dataNormt <- t(dataNorm)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="euclidean")
#groups <- factor(c(rep(1,16), rep(2,8))
groups <- c(metaD$Group)
#mod <- betadisper(dist.Metab ~ Condition + Time, metaD)
mod <- betadisper(dist.Metab, groups)
#Calculates spatial median, not center of mass
#plot(mod)
#pcoa.Metab <- cmdscale(dist.Metab, eig=TRUE)
#Variance explained should match individual plots from Metaboanalyst
PC1var <- mod$eig[1]/sum(mod$eig) #Checks out
PC2var <- mod$eig[2]/sum(mod$eig) #Checks out

#Extract scores

modScores <- scores(mod)
#Extract centroid
centroids <- as.data.frame(modScores$centroids)
#Extract species and time info
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
#Remove redundancies
condI <- condI[!duplicated(condI),]
#Add species
centroids$Label <- condI$Condition
#Add time
centroids$Time <- condI$Time
#Flip y-axis
centroids$PCoA2 <- centroids$PCoA2*-1


#Extract axes vectors from individual samples
sites <- as.data.frame(modScores$sites)
#add groups
sites$groups <- mod$group


#Calculated PCoA axes SDs for each group
ag <- aggregate(.~groups,sites, function(x) sd=sd(x))
#Std axis 1
sd_axis1 <- ag$PCoA1
#Std axis 2
sd_axis2 <- ag$PCoA2

library(ggplot2)
#library(repr)
#options(repr.plot.width = 5, repr.plot.height = 3)
#Change time from factor to numeric
centroids$Time <- rep(c(12.5,25,30,35,40,45),3)

PCA_PolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% var. explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% var. explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 20),legend.title=element_text(size=20))

PCA_PN_Bt <- PCA_PolarNeg + ylim(-125,0) + xlim(-250,-125) + theme(axis.title=element_blank())
PCA_PN_Cv <- PCA_PolarNeg + ylim(80,190) + xlim(15,60) + theme(axis.title=element_blank())
PCA_PN_Ps <- PCA_PolarNeg + ylim(-150,-50) + xlim(75,180) + theme(axis.title=element_blank())

library(patchwork)
#PCA_PolarNeg_plots <- wrap_plots(PCA_PolarNeg,PCA_PN_Cv,PCA_PN_Bt,PCA_PN_Ps)
#  + theme(legend.position="right") + guides(fill=guide_legend(title = "Strain",order=1),size=guide_legend(title = "Time",order=2),color=guide_legend(order=3))

PCA_PolarNeg_plots <- PCA_PolarNeg + PCA_PN_Cv +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_PN_Bt + PCA_PN_Ps + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')

ggsave("Figures/Fig1/PCA_PolarNeg.eps",plot=PCA_PolarNeg_plots,device="eps",width=15,height=10,dpi=600)

#Run variation partitioning.
PNvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.66053 #Time: 0.00160 #Species + Time:0.71284





#####NonPolar Pos#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Species"=mDataSplit$X2, "Time"=mDataSplitTime$X1)

#Prepare metadata table for varpart
metaDvP <- metaD
#Add row names
row.names(metaDvP) <- mData

metaD <- data.frame("Group"=paste(mDataSplit$X2,mDataSplitTime$X1))
#Remove column names, remove row 1, and transpose dataframe
#colnames(dataNorm) <- NULL
dataNorm <- dataNorm[-1,]
#Convert characters to numeric
dataNorm <- as.data.frame(lapply(dataNorm, as.numeric))

#transpose
dataNormt <- t(dataNorm)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="euclidean")
#groups <- factor(c(rep(1,16), rep(2,8))
groups <- c(metaD$Group)
#mod <- betadisper(dist.Metab ~ Condition + Time, metaD)
mod <- betadisper(dist.Metab, groups)
#Calculates spatial median, not center of mass
#plot(mod)
#pcoa.Metab <- cmdscale(dist.Metab, eig=TRUE)
#Variance explained should match individual plots from Metaboanalyst
PC1var <- mod$eig[1]/sum(mod$eig) #Checks out
PC2var <- mod$eig[2]/sum(mod$eig) #Checks out

#Extract scores

modScores <- scores(mod)
#Extract centroid
centroids <- as.data.frame(modScores$centroids)
#Extract species and time info
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
#Remove redundancies
condI <- condI[!duplicated(condI),]
#Add species
centroids$Label <- condI$Condition
#Add time
centroids$Time <- condI$Time
#Flip x-axis
centroids$PCoA1 <- centroids$PCoA1*-1


#Extract axes vectors from individual samples
sites <- as.data.frame(modScores$sites)
#add groups
sites$groups <- mod$group


#Calculated PCoA axes SDs for each group
ag <- aggregate(.~groups,sites, function(x) sd=sd(x))
#Std axis 1
sd_axis1 <- ag$PCoA1
#Std axis 2
sd_axis2 <- ag$PCoA2

library(ggplot2)
#library(repr)
#options(repr.plot.width = 5, repr.plot.height = 3)
#Change time from factor to numeric
centroids$Time <- rep(c(12.5,25,30,35,40,45),3)

PCA_NonPolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% var. explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% var. explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position="none",axis.title = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_text(size=20))

PCA_NPP_Bt <- PCA_NonPolarPos + ylim(-70,-37.5) + xlim(-40,-17.5) + theme(axis.title=element_blank())
PCA_NPP_Cv <- PCA_NonPolarPos + ylim(30,70) + xlim(-45,-10) + theme(axis.title=element_blank())
PCA_NPP_Ps <- PCA_NonPolarPos + ylim(-20,25) + xlim(45,110) + theme(axis.title=element_blank())

library(patchwork)
#PCA_NonPolarPos_plots <- wrap_plots(PCA_NonPolarPos,PCA_NPP_Cv,PCA_NPP_Bt,PCA_NPP_Ps)
#ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos_plots,device="png",width=15,height=8.7,dpi=300)

PCA_NonPolarPos_plots <- PCA_NonPolarPos + PCA_NPP_Cv +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_NPP_Bt + PCA_NPP_Ps + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')

ggsave("Figures/Fig1/PCA_NonPolarPos.eps",plot=PCA_NonPolarPos_plots,device="eps",width=15,height=10,dpi=600)


#Run variation partitioning.
NPPvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.80584 #Time: -0.03962 #Species + Time:0.83962





#####NonPolar Neg#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Species"=mDataSplit$X2, "Time"=mDataSplitTime$X1)

#Prepare metadata table for varpart
metaDvP <- metaD
#Add row names
row.names(metaDvP) <- mData

metaD <- data.frame("Group"=paste(mDataSplit$X2,mDataSplitTime$X1))
#Remove column names, remove row 1, and transpose dataframe
#colnames(dataNorm) <- NULL
dataNorm <- dataNorm[-1,]
#Convert characters to numeric
dataNorm <- as.data.frame(lapply(dataNorm, as.numeric))

#transpose
dataNormt <- t(dataNorm)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="euclidean")
#groups <- factor(c(rep(1,16), rep(2,8))
groups <- c(metaD$Group)
#mod <- betadisper(dist.Metab ~ Condition + Time, metaD)
mod <- betadisper(dist.Metab, groups)
#Calculates spatial median, not center of mass
#plot(mod)
#pcoa.Metab <- cmdscale(dist.Metab, eig=TRUE)
#Variance explained should match individual plots from Metaboanalyst
PC1var <- mod$eig[1]/sum(mod$eig) #Checks out
PC2var <- mod$eig[2]/sum(mod$eig) #Checks out

#Extract scores

modScores <- scores(mod)
#Extract centroid
centroids <- as.data.frame(modScores$centroids)
#Extract species and time info
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
#Remove redundancies
condI <- condI[!duplicated(condI),]
#Add species
centroids$Label <- condI$Condition
#Add time
centroids$Time <- condI$Time
#Flip y-axis
centroids$PCoA2 <- centroids$PCoA2*-1


#Extract axes vectors from individual samples
sites <- as.data.frame(modScores$sites)
#add groups
sites$groups <- mod$group


#Calculated PCoA axes SDs for each group
ag <- aggregate(.~groups,sites, function(x) sd=sd(x))
#Std axis 1
sd_axis1 <- ag$PCoA1
#Std axis 2
sd_axis2 <- ag$PCoA2

library(ggplot2)
#library(repr)
#options(repr.plot.width = 5, repr.plot.height = 3)
#Change time from factor to numeric
centroids$Time <- rep(c(12.5,25,30,35,40,45),3)

PCA_NonPolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% var. explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% var. explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 20),legend.title=element_text(size=20))

PCA_NPN_Bt <- PCA_NonPolarNeg + ylim(5,27.5) + xlim(-100,-55) + theme(axis.title=element_blank())
PCA_NPN_Cv <- PCA_NonPolarNeg + ylim(22.5,47.5) + xlim(45,65) + theme(axis.title=element_blank())
PCA_NPN_Ps <- PCA_NonPolarNeg + ylim(-90,-27.5) + xlim(-15,55) + theme(axis.title=element_blank())

library(patchwork)
#PCA_NonPolarNeg_plots <- wrap_plots(PCA_NonPolarNeg,PCA_NPN_Cv,PCA_NPN_Bt,PCA_NPN_Ps)

PCA_NonPolarNeg_plots <- PCA_NonPolarNeg + PCA_NPN_Cv +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_NPN_Bt + PCA_NPN_Ps + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')

ggsave("Figures/Fig1/PCA_NonPolarNeg.eps",plot=PCA_NonPolarNeg_plots,device="eps",width=15,height=10,dpi=600)


#Run variation partitioning.
NPNvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.81520 #Time: -0.04966 #Species + Time:0.85304



#Let's put all mass spec analyses together

library(patchwork)

PCA_plots <- PCA_PolarPos + PCA_PolarNeg +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_NonPolarPos + PCA_NonPolarNeg + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')


ggsave("Figures/Fig1/PCA_Plots.eps",plot=PCA_plots,device="eps",width=15,height=10,dpi=600)


#Perform Protest
library(vegan)
library(dplyr)

#####Polar Pos#####

data <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarPos_pcaScores.csv", header=TRUE)

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(data, Strain == "Bt" & Name == "R1")
BtR2 <- filter(data, Strain == "Bt" & Name == "R2")
BtR3 <- filter(data, Strain == "Bt" & Name == "R3")
BtR4 <- filter(data, Strain == "Bt" & Name == "R4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(data, Strain == "Cv" & Name == "R1")
CvR2 <- filter(data, Strain == "Cv" & Name == "R2")
CvR3 <- filter(data, Strain == "Cv" & Name == "R3")
CvR4 <- filter(data, Strain == "Cv" & Name == "R4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(data, Strain == "Ps" & Name == "R1")
PsR2 <- filter(data, Strain == "Ps" & Name == "R2")
PsR3 <- filter(data, Strain == "Ps" & Name == "R3")
PsR4 <- filter(data, Strain == "Ps" & Name == "R4")

#Perform protest for Bt
protest(X=BtR1[,2:3],Y=BtR2[,2:3])
protest(X=BtR1[,2:3],Y=BtR3[,2:3])
protest(X=BtR1[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR2[,2:3],Y=BtR3[,2:3])
protest(X=BtR2[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,2:3],Y=CvR2[,2:3])
protest(X=CvR1[,2:3],Y=CvR3[,2:3])
protest(X=CvR1[,2:3],Y=CvR4[,2:3])
protest(X=CvR2[,2:3],Y=CvR3[,2:3])
protest(X=CvR2[,2:3],Y=CvR4[,2:3])
protest(X=CvR3[,2:3],Y=CvR4[,2:3])

#Perform protest for Ps
protest(X=PsR1[,2:3],Y=PsR2[,2:3])
protest(X=PsR1[,2:3],Y=PsR3[,2:3])
protest(X=PsR1[,2:3],Y=PsR4[,2:3])
protest(X=PsR2[,2:3],Y=PsR3[,2:3])
protest(X=PsR2[,2:3],Y=PsR4[,2:3])
protest(X=PsR3[,2:3],Y=PsR4[,2:3])

#####Polar Neg#####

data <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarNeg_pcaScores.csv", header=TRUE)

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(data, Strain == "Bt" & Name == "R1")
BtR2 <- filter(data, Strain == "Bt" & Name == "R2")
BtR3 <- filter(data, Strain == "Bt" & Name == "R3")
BtR4 <- filter(data, Strain == "Bt" & Name == "R4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(data, Strain == "Cv" & Name == "R1")
CvR2 <- filter(data, Strain == "Cv" & Name == "R2")
CvR3 <- filter(data, Strain == "Cv" & Name == "R3")
CvR4 <- filter(data, Strain == "Cv" & Name == "R4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(data, Strain == "Ps" & Name == "R1")
PsR2 <- filter(data, Strain == "Ps" & Name == "R2")
PsR3 <- filter(data, Strain == "Ps" & Name == "R3")
PsR4 <- filter(data, Strain == "Ps" & Name == "R4")

#Perform protest for Bt
protest(X=BtR1[,2:3],Y=BtR2[,2:3])
protest(X=BtR1[,2:3],Y=BtR3[,2:3])
protest(X=BtR1[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR2[,2:3],Y=BtR3[,2:3])
protest(X=BtR2[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,2:3],Y=CvR2[,2:3])
protest(X=CvR1[,2:3],Y=CvR3[,2:3])
protest(X=CvR1[,2:3],Y=CvR4[,2:3])
protest(X=CvR2[,2:3],Y=CvR3[,2:3])
protest(X=CvR2[,2:3],Y=CvR4[,2:3])
protest(X=CvR3[,2:3],Y=CvR4[,2:3])

#Perform protest for Ps
protest(X=PsR1[,2:3],Y=PsR2[,2:3])
protest(X=PsR1[,2:3],Y=PsR3[,2:3])
protest(X=PsR1[,2:3],Y=PsR4[,2:3])
protest(X=PsR2[,2:3],Y=PsR3[,2:3])
protest(X=PsR2[,2:3],Y=PsR4[,2:3])
protest(X=PsR3[,2:3],Y=PsR4[,2:3])

#####NonPolar Pos#####

data <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarPos_pcaScores.csv", header=TRUE)

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(data, Strain == "Bt" & Name == "R1")
BtR2 <- filter(data, Strain == "Bt" & Name == "R2")
BtR3 <- filter(data, Strain == "Bt" & Name == "R3")
BtR4 <- filter(data, Strain == "Bt" & Name == "R4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(data, Strain == "Cv" & Name == "R1")
CvR2 <- filter(data, Strain == "Cv" & Name == "R2")
CvR3 <- filter(data, Strain == "Cv" & Name == "R3")
CvR4 <- filter(data, Strain == "Cv" & Name == "R4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(data, Strain == "Ps" & Name == "R1")
PsR2 <- filter(data, Strain == "Ps" & Name == "R2")
PsR3 <- filter(data, Strain == "Ps" & Name == "R3")
PsR4 <- filter(data, Strain == "Ps" & Name == "R4")

#Perform protest for Bt- not enough data for R1
protest(X=BtR2[,2:3],Y=BtR3[,2:3])
protest(X=BtR2[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,2:3],Y=CvR2[2:6,2:3]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR1[,2:3],Y=CvR3[2:6,2:3]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR1[,2:3],Y=CvR4[2:6,2:3]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR2[,2:3],Y=CvR3[,2:3])
protest(X=CvR2[,2:3],Y=CvR4[,2:3])
protest(X=CvR3[,2:3],Y=CvR4[,2:3])

#Perform protest for Ps - no data for R1
protest(X=PsR2[,2:3],Y=PsR3[,2:3])
protest(X=PsR2[2:6,2:3],Y=PsR4[,2:3]) #Remove first time point- missing data for Ps12.5h_BioRep4
protest(X=PsR3[2:6,2:3],Y=PsR4[,2:3]) #Remove first time point- missing data for Ps12.5h_BioRep4

#####NonPolar Neg#####

data <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarNeg_pcaScores.csv", header=TRUE)

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(data, Strain == "Bt" & Name == "R1")
BtR2 <- filter(data, Strain == "Bt" & Name == "R2")
BtR3 <- filter(data, Strain == "Bt" & Name == "R3")
BtR4 <- filter(data, Strain == "Bt" & Name == "R4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(data, Strain == "Cv" & Name == "R1")
CvR2 <- filter(data, Strain == "Cv" & Name == "R2")
CvR3 <- filter(data, Strain == "Cv" & Name == "R3")
CvR4 <- filter(data, Strain == "Cv" & Name == "R4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(data, Strain == "Ps" & Name == "R1")
PsR2 <- filter(data, Strain == "Ps" & Name == "R2")
PsR3 <- filter(data, Strain == "Ps" & Name == "R3")
PsR4 <- filter(data, Strain == "Ps" & Name == "R4")

#Perform protest for Bt- not enough data for R1
protest(X=BtR2[1:5,2:3],Y=BtR3[,2:3]) #Remove last time point- missing data for Bt45h_BioRep3
protest(X=BtR2[1:5,2:3],Y=BtR4[,2:3]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[,2:3],Y=BtR4[,2:3])

#Perform protest for Cv- not enough data for R1
protest(X=CvR2[,2:3],Y=CvR3[,2:3])
protest(X=CvR2[,2:3],Y=CvR4[,2:3])
protest(X=CvR3[,2:3],Y=CvR4[,2:3])

#Perform protest for Ps - no data for R1
protest(X=PsR2[,2:3],Y=PsR3[,2:3])
protest(X=PsR2[2:6,2:3],Y=PsR4[,2:3]) #Remove first time point- missing data for Ps12.5h_BioRep4
protest(X=PsR3[2:6,2:3],Y=PsR4[,2:3]) #Remove first time point- missing data for Ps12.5h_BioRep4






###summary
all r2 > 0.9383
all p <= 0.025





#Make panels sizes the same size and save plots
#library(gridExtra)
#library(egg)

#pca_fixedPP <- set_panel_size(PCA_PolarPos,  width  = unit(5, "in"), height = unit(3, "in"))
#pca_fixedPN <- set_panel_size(PCA_PolarNeg,  width  = unit(5, "in"), height = unit(3, "in"))
#pca_fixedNPP <- set_panel_size(PCA_NonPolarPos,  width  = unit(15, "in"), height = unit(8.7, "in"))
#pca_fixedNPN <- set_panel_size(PCA_NonPolarNeg,  width  = unit(15, "in"), height = unit(8.7, "in"))

#gaPP <- arrangeGrob(pca_fixedPP)
#gaPN <- arrangeGrob(pca_fixedPN)
#gaNPP <- arrangeGrob(pca_fixedNPP)
#gaNPN <- arrangeGrob(pca_fixedNPN)

#ggsave("Figures/Fig1/PCA_PolarPos.png",plot=gaPP,device="png",width=7,height=5,dpi=300)
#ggsave("Figures/Fig1/PCA_PolarNeg.png",plot=gaPN,device="png",width=7,height=5,dpi=300)
#ggsave("Figures/Fig1/PCA_NonPolarPos.png",plot=gaNPP,device="png",width=20,height=11.6,dpi=600)
#ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=gaNPN,device="png",width=20,height=11.6,dpi=600)
