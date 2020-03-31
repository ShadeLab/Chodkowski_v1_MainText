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
PCA_PolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "none")

PCA_PP_Bt <- PCA_PolarPos + ylim(-55,25) + xlim(-200,-80) + theme(axis.title=element_blank())
PCA_PP_Cv <- PCA_PolarPos + ylim(47.5,60) + xlim(42.5,57) + theme(axis.title=element_blank())
PCA_PP_Ps <- PCA_PolarPos + ylim(-80,-30) + xlim(50,90) + theme(axis.title=element_blank())

library(patchwork)
PCA_PolarPos_plots <- wrap_plots(PCA_PolarPos,PCA_PP_Cv,PCA_PP_Bt,PCA_PP_Ps)
ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos_plots,device="png",width=15,height=8.7,dpi=300)

  #theme_bw() +
  #theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos,device="png",width=15,height=8.7,dpi=600)

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
PCA_PolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "none")

PCA_PN_Bt <- PCA_PolarNeg + ylim(-125,0) + xlim(-250,-125) + theme(axis.title=element_blank())
PCA_PN_Cv <- PCA_PolarNeg + ylim(80,190) + xlim(15,60) + theme(axis.title=element_blank())
PCA_PN_Ps <- PCA_PolarNeg + ylim(-150,-50) + xlim(75,180) + theme(axis.title=element_blank())

library(patchwork)
PCA_PolarNeg_plots <- wrap_plots(PCA_PolarNeg,PCA_PN_Cv,PCA_PN_Bt,PCA_PN_Ps)
ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos_plots,device="png",width=15,height=8.7,dpi=300)

  #theme_bw() +
  #theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_PolarNeg.png",plot=PCA_PolarNeg,device="png",width=15,height=8.7,units="in",dpi=600)

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
PCA_NonPolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "none")

PCA_PP_Bt <- PCA_PolarPos + ylim(-55,25) + xlim(-200,-80) + theme(axis.title=element_blank())
PCA_PP_Cv <- PCA_PolarPos + ylim(47.5,60) + xlim(42.5,57) + theme(axis.title=element_blank())
PCA_PP_Ps <- PCA_PolarPos + ylim(-80,-30) + xlim(50,90) + theme(axis.title=element_blank())

library(patchwork)
PCA_PolarPos_plots <- wrap_plots(PCA_PolarPos,PCA_PP_Cv,PCA_PP_Bt,PCA_PP_Ps)
ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos_plots,device="png",width=15,height=8.7,dpi=300)

  #theme_bw() +
  #theme(legend.position = "right")


#ggsave("Figures/Fig1/PCA_NonPolarPos.png",plot=PCA_NonPolarPos,device="png",width=15,height=8.7,dpi=600)

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
PCA_NonPolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  #coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(2,3.5,5,7,9,12))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "none")

PCA_NPN_Bt <- PCA_NonPolarNeg + ylim(-125,0) + xlim(-250,-125) + theme(axis.title=element_blank())
PCA_NPN_Cv <- PCA_NonPolarNeg + ylim(80,190) + xlim(15,60) + theme(axis.title=element_blank())
PCA_NPN_Ps <- PCA_NonPolarNeg + ylim(-150,-50) + xlim(75,180) + theme(axis.title=element_blank())

library(patchwork)
PCA_NonPolarNeg_plots <- wrap_plots(PCA_NonPolarNeg,PCA_NPN_Cv,PCA_NPN_Bt,PCA_NPN_Ps)
ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=PCA_NonPolarNeg_plots,device="png",width=15,height=8.7,dpi=300)

  #theme_bw() +
  #theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=PCA_NonPolarNeg,device="png",width=15,height=8.7,dpi=600)

#Run variation partitioning.
NPNvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.81520 #Time: -0.04966 #Species + Time:0.85304



#Let's put all mass spec analyses together

#Make panels sizes the same size and save plots
library(gridExtra)
library(egg)

pca_fixedPP <- set_panel_size(PCA_PolarPos,  width  = unit(5, "in"), height = unit(3, "in"))
pca_fixedPN <- set_panel_size(PCA_PolarNeg,  width  = unit(5, "in"), height = unit(3, "in"))
pca_fixedNPP <- set_panel_size(PCA_NonPolarPos,  width  = unit(15, "in"), height = unit(8.7, "in"))
pca_fixedNPN <- set_panel_size(PCA_NonPolarNeg,  width  = unit(15, "in"), height = unit(8.7, "in"))

gaPP <- arrangeGrob(pca_fixedPP)
gaPN <- arrangeGrob(pca_fixedPN)
gaNPP <- arrangeGrob(pca_fixedNPP)
gaNPN <- arrangeGrob(pca_fixedNPN)

ggsave("Figures/Fig1/PCA_PolarPos.png",plot=gaPP,device="png",width=7,height=5,dpi=300)
ggsave("Figures/Fig1/PCA_PolarNeg.png",plot=gaPN,device="png",width=7,height=5,dpi=300)
ggsave("Figures/Fig1/PCA_NonPolarPos.png",plot=gaNPP,device="png",width=20,height=11.6,dpi=600)
ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=gaNPN,device="png",width=20,height=11.6,dpi=600)
