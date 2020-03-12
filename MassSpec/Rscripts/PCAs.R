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
metaD <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
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
plot(mod)
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

PCA_PolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  theme_bw() +
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(0.1,0.75,1.5,3,4,5))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_PolarPos.png",plot=PCA_PolarPos,device="png",width=15,height=8.7,dpi=600)

#####Polar Neg#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
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
plot(mod)
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

PCA_PolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  theme_bw() +
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(0.1,0.75,1.5,3,4,5))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_PolarNeg.png",plot=PCA_PolarNeg,device="png",width=15,height=8.7,units="in",dpi=600)

#####NonPolar Pos#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
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
plot(mod)
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

PCA_NonPolarPos <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  theme_bw() +
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(0.1,0.75,1.5,3,4,5))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_NonPolarPos.png",plot=PCA_NonPolarPos,device="png",width=15,height=8.7,dpi=600)

#####NonPolar Neg#####
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
#Create metaData from column headers
mData <- colnames(dataNorm)

#Split names by "_" to get species
mDataSplit <- data.frame(do.call('rbind', strsplit(as.character(mData),'_',fixed=TRUE)))
#Split time remove replicate
mDataSplitTime <- data.frame(do.call('rbind', strsplit(as.character(mDataSplit$X3),'.',fixed=TRUE)))
#Bring metaData together
metaD <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
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
plot(mod)
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

PCA_NonPolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill= factor(Label),size=factor(Time)), colour="black",shape=21)+
  coord_fixed()+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1+sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1-sd_axis1,y=PCoA2,yend=PCoA2,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2+sd_axis1,color=Label))+
  geom_segment(data=centroids, aes(x=PCoA1,xend=PCoA1,y=PCoA2,yend=PCoA2-sd_axis1,color=Label))+
  xlab(label = paste("PC1"," (", round(PC1var,digits = 3)*100, "% variance explained)", sep = ""))+
  ylab(label = paste("PC2"," (", round(PC2var,digits = 3)*100, "% variance explained)", sep = ""))+
  theme_bw() +
  scale_fill_manual(values=c("#56B4E9", "#9900CC","#33CC00")) +
  scale_size_manual(values = c(0.1,0.75,1.5,3,4,5))+
  scale_color_manual(values=c("#56B4E9", "#9900CC","#33CC00"))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 7)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  theme(legend.position = "right")

#ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=PCA_NonPolarNeg,device="png",width=15,height=8.7,dpi=600)


#Make panels sizes the same size and save plots
library(gridExtra)
library(egg)

pca_fixedPP <- set_panel_size(PCA_PolarPos,  width  = unit(15, "in"), height = unit(8.7, "in"))
pca_fixedPN <- set_panel_size(PCA_PolarNeg,  width  = unit(15, "in"), height = unit(8.7, "in"))
pca_fixedNPP <- set_panel_size(PCA_NonPolarPos,  width  = unit(15, "in"), height = unit(8.7, "in"))
pca_fixedNPN <- set_panel_size(PCA_NonPolarNeg,  width  = unit(15, "in"), height = unit(8.7, "in"))

gaPP <- arrangeGrob(pca_fixedPP)
gaPN <- arrangeGrob(pca_fixedPN)
gaNPP <- arrangeGrob(pca_fixedNPP)
gaNPN <- arrangeGrob(pca_fixedNPN)

ggsave("Figures/Fig1/PCA_PolarPos.png",plot=gaPP,device="png",width=20,height=11.6,dpi=600)
ggsave("Figures/Fig1/PCA_PolarNeg.png",plot=gaPN,device="png",width=20,height=11.6,dpi=600)
ggsave("Figures/Fig1/PCA_NonPolarPos.png",plot=gaNPP,device="png",width=20,height=11.6,dpi=600)
ggsave("Figures/Fig1/PCA_NonPolarNeg.png",plot=gaNPN,device="png",width=20,height=11.6,dpi=600)
