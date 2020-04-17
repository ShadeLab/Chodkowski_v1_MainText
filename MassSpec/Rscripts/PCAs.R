#Read in data

###Polar Pos###

#Read in NA indicies
NA_Index <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PolarPos_NA_Index.rds")

#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Convert indicies of interest to NAs. These features did not pass release criteria for particular strains.
for(i in 1:nrow(NA_Index)){
dataNorm.edit[NA_Index$row[i],NA_Index$col[i]] <- NA
}

#Create metaData from column headers
mData <- colnames(dataNorm.edit)

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

#Convert characters to numeric
dataNorm.edit <- as.data.frame(lapply(dataNorm.edit, as.numeric))

#Transpose
dataNormt <- t(dataNorm.edit)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="bray")

#Creat groups
groups <- c(metaD$Group)
#Calculate betadisper
mod <- betadisper(dist.Metab, groups)
#Calculates spatial median, not center of mass

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
#centroids$PCoA2 <- centroids$PCoA2*-1

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
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 5)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 5)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title=element_text(size=20))

#Run variation partitioning.
PPvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#####Varpart when non-released replaced by NAs
#Species: 0.35787 #Time: 0.23528 #Species x Time:0.62992
#####Varpart when non-released replaced by 0s
#Species: 0.89178 #Time: -0.03724 #Species x Time:0.92299

#Run adonis
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
rownames(condI) <- rownames(dataNormt)
adonis(dist.Metab ~ Condition*Time, data=condI,permutations=999)

#Perform pairwise adonis post-hoc
library(pairwiseAdonis)
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.adonis(dist.Metab,condI$Cond)

#####Perform Protest#####
library(vegan)
library(dplyr)

#Extract PC1&2 for bio reps
mod_sample.PCs <- as.data.frame(modScores$sites)
#Add strains 
mod_sample.PCs$Strain <- mDataSplit$X2
#Add bioreps
mod_sample.PCs$BR <- mDataSplitTime$X2

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "1")
BtR2 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "2")
BtR3 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "3")
BtR4 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "1")
CvR2 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "2")
CvR3 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "3")
CvR4 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "1")
PsR2 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "2")
PsR3 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "3")
PsR4 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "4")

#Perform protest for Bt
protest(X=BtR1[,1:2],Y=BtR2[,1:2])
protest(X=BtR1[,1:2],Y=BtR3[,1:2])
protest(X=BtR1[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR2[,1:2],Y=BtR3[,1:2])
protest(X=BtR2[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,1:2],Y=CvR2[,1:2])
protest(X=CvR1[,1:2],Y=CvR3[,1:2])
protest(X=CvR1[,1:2],Y=CvR4[,1:2])
protest(X=CvR2[,1:2],Y=CvR3[,1:2])
protest(X=CvR2[,1:2],Y=CvR4[,1:2])
protest(X=CvR3[,1:2],Y=CvR4[,1:2])

#Perform protest for Ps
protest(X=PsR1[,1:2],Y=PsR2[,1:2])
protest(X=PsR1[,1:2],Y=PsR3[,1:2])
protest(X=PsR1[,1:2],Y=PsR4[,1:2])
protest(X=PsR2[,1:2],Y=PsR3[,1:2])
protest(X=PsR2[,1:2],Y=PsR4[,1:2])
protest(X=PsR3[,1:2],Y=PsR4[,1:2])



#####Polar Neg#####
#Read in NA indicies
NA_Index <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PolarNeg_NA_Index.rds")

dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Convert indicies of interest to NAs. These features did not pass release criteria for particular strains.
for(i in 1:nrow(NA_Index)){
dataNorm.edit[NA_Index$row[i],NA_Index$col[i]] <- NA
}

#Create metaData from column headers
mData <- colnames(dataNorm.edit)

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

#Convert characters to numeric
dataNorm.edit <- as.data.frame(lapply(dataNorm.edit, as.numeric))

#transpose
dataNormt <- t(dataNorm.edit)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="bray")
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
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 5)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 5)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title=element_text(size=20))

#Run variation partitioning.
PNvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.37597 #Time: 0.18233 #Species x Time:0.59106
#####Varpart when non-released replaced by 0s
#Species: 0.88964 #Time: -0.04022 #Species x Time:0.91744

#Run adonis
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
rownames(condI) <- rownames(dataNormt)
adonis(dist.Metab ~ Condition*Time, data=condI,permutations=999)

#Perform pairwise adonis post-hoc
library(pairwiseAdonis)
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.adonis(dist.Metab,condI$Cond)

#####Perform Protest#####
library(vegan)
library(dplyr)

#Extract PC1&2 for bio reps
mod_sample.PCs <- as.data.frame(modScores$sites)
#Add strains 
mod_sample.PCs$Strain <- mDataSplit$X2
#Add bioreps
mod_sample.PCs$BR <- mDataSplitTime$X2

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "1")
BtR2 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "2")
BtR3 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "3")
BtR4 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "1")
CvR2 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "2")
CvR3 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "3")
CvR4 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "1")
PsR2 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "2")
PsR3 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "3")
PsR4 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "4")

#Perform protest for Bt
protest(X=BtR1[,1:2],Y=BtR2[,1:2])
protest(X=BtR1[,1:2],Y=BtR3[,1:2])
protest(X=BtR1[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR2[,1:2],Y=BtR3[,1:2])
protest(X=BtR2[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,1:2],Y=CvR2[,1:2])
protest(X=CvR1[,1:2],Y=CvR3[,1:2])
protest(X=CvR1[,1:2],Y=CvR4[,1:2])
protest(X=CvR2[,1:2],Y=CvR3[,1:2])
protest(X=CvR2[,1:2],Y=CvR4[,1:2])
protest(X=CvR3[,1:2],Y=CvR4[,1:2])

#Perform protest for Ps
protest(X=PsR1[,1:2],Y=PsR2[,1:2])
protest(X=PsR1[,1:2],Y=PsR3[,1:2])
protest(X=PsR1[,1:2],Y=PsR4[,1:2])
protest(X=PsR2[,1:2],Y=PsR3[,1:2])
protest(X=PsR2[,1:2],Y=PsR4[,1:2])
protest(X=PsR3[,1:2],Y=PsR4[,1:2])




#####NonPolar Pos#####
NA_Index <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/NonPolarPos_NA_Index.rds")

#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Convert indicies of interest to NAs. These features did not pass release criteria for particular strains.
for(i in 1:nrow(NA_Index)){
dataNorm.edit[NA_Index$row[i],NA_Index$col[i]] <- NA
}

#Create metaData from column headers
mData <- colnames(dataNorm.edit)

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

#Convert characters to numeric
dataNorm.edit <- as.data.frame(lapply(dataNorm.edit, as.numeric))

#transpose
dataNormt <- t(dataNorm.edit)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="bray")
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
#centroids$PCoA1 <- centroids$PCoA1*-1


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
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 5)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 5)) +
  theme(legend.position="none",axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.title=element_text(size=20))

#Run variation partitioning.
NPPvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.31800 #Time: 0.38165 #Species x Time:0.75185
#####Varpart when non-released replaced by 0s
#Species: 0.91854 #Time: -0.05818 #Species x Time:0.94454

#Run adonis
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
rownames(condI) <- rownames(dataNormt)
adonis(dist.Metab ~ Condition*Time, data=condI,permutations=999)

#Perform pairwise adonis post-hoc
library(pairwiseAdonis)
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.adonis(dist.Metab,condI$Cond)

#####Perform Protest#####
library(vegan)
library(dplyr)

#Extract PC1&2 for bio reps
mod_sample.PCs <- as.data.frame(modScores$sites)
#Add strains 
mod_sample.PCs$Strain <- mDataSplit$X2
#Add bioreps
mod_sample.PCs$BR <- mDataSplitTime$X2

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "1")
BtR2 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "2")
BtR3 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "3")
BtR4 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "1")
CvR2 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "2")
CvR3 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "3")
CvR4 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "1")
PsR2 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "2")
PsR3 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "3")
PsR4 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "4")

#Perform protest for Bt- not enough data for R1
protest(X=BtR2[,1:2],Y=BtR3[,1:2])
protest(X=BtR2[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4

#Perform protest for Cv
protest(X=CvR1[,1:2],Y=CvR2[2:6,1:2]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR1[,1:2],Y=CvR3[2:6,1:2]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR1[,1:2],Y=CvR4[2:6,1:2]) #Remove first time point- missing data for Cv12.5h_BioRep1
protest(X=CvR2[,1:2],Y=CvR3[,1:2])
protest(X=CvR2[,1:2],Y=CvR4[,1:2])
protest(X=CvR3[,1:2],Y=CvR4[,1:2])

#Perform protest for Ps - no data for R1
protest(X=PsR2[,1:2],Y=PsR3[,1:2])
protest(X=PsR2[2:6,1:2],Y=PsR4[,1:2]) #Remove first time point- missing data for Ps12.5h_BioRep4
protest(X=PsR3[2:6,1:2],Y=PsR4[,1:2]) #Remove first time point- missing data for Ps12.5h_BioRep4




#####NonPolar Neg#####
NA_Index <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/NonPolarNeg_NA_Index.rds")

#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Convert indicies of interest to NAs. These features did not pass release criteria for particular strains.
for(i in 1:nrow(NA_Index)){
dataNorm.edit[NA_Index$row[i],NA_Index$col[i]] <- NA
}

#Create metaData from column headers
mData <- colnames(dataNorm.edit)

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

#Convert characters to numeric
dataNorm.edit <- as.data.frame(lapply(dataNorm.edit, as.numeric))

#transpose
dataNormt <- t(dataNorm.edit)

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="bray")
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
#centroids$PCoA2 <- centroids$PCoA2*-1


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
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 5)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 5)) +
  theme(legend.position="none",axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title=element_text(size=20))

#Run variation partitioning.
NPNvarPart <- varpart(dist.Metab,~Species,~Time,data=metaDvP)
#Species: 0.34374 #Time: 0.40226 #Species x Time:0.79137
#####Varpart when non-released replaced by 0s
#Species: 0.91629 #Time: -0.06694 #Species x Time:0.94001

#Run adonis
condI <- data.frame("Condition"=mDataSplit$X2, "Time"=mDataSplitTime$X1)
rownames(condI) <- rownames(dataNormt)
adonis(dist.Metab ~ Condition*Time, data=condI,permutations=999)

#Perform pairwise adonis post-hoc
library(pairwiseAdonis)
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.adonis(dist.Metab,condI$Cond)

#####Perform Protest#####
library(vegan)
library(dplyr)

#Extract PC1&2 for bio reps
mod_sample.PCs <- as.data.frame(modScores$sites)
#Add strains 
mod_sample.PCs$Strain <- mDataSplit$X2
#Add bioreps
mod_sample.PCs$BR <- mDataSplitTime$X2

#Extract B. thailandensis PCA PC1 & PC2 scores
BtR1 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "1")
BtR2 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "2")
BtR3 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "3")
BtR4 <- filter(mod_sample.PCs, Strain == "Bt" & BR == "4")

#Extract C. violaceum PCA PC1 & PC2 scores
CvR1 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "1")
CvR2 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "2")
CvR3 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "3")
CvR4 <- filter(mod_sample.PCs, Strain == "Cv" & BR == "4")

#Extract P. syringae PCA PC1 & PC2 scores
PsR1 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "1")
PsR2 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "2")
PsR3 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "3")
PsR4 <- filter(mod_sample.PCs, Strain == "Ps" & BR == "4")

#Perform protest for Bt- not enough data for R1
protest(X=BtR2[,1:2],Y=BtR3[,1:2]) #Check if pseudo sample was added
protest(X=BtR2[1:5,1:2],Y=BtR4[,1:2]) #Remove last time point- missing data for Bt45h_BioRep4
protest(X=BtR3[1:5,1:2],Y=BtR4[,1:2])

#Perform protest for Cv- not enough data for R1
protest(X=CvR2[,1:2],Y=CvR3[,1:2])
protest(X=CvR2[,1:2],Y=CvR4[,1:2])
protest(X=CvR3[,1:2],Y=CvR4[,1:2])

#Perform protest for Ps - no data for R1
protest(X=PsR2[,1:2],Y=PsR3[,1:2])
protest(X=PsR2[2:6,1:2],Y=PsR4[,1:2]) #Remove first time point- missing data for Ps12.5h_BioRep4
protest(X=PsR3[2:6,1:2],Y=PsR4[,1:2]) #Remove first time point- missing data for Ps12.5h_BioRep4




#Let's put all mass spec analyses together

library(patchwork)

PCA_plots <- PCA_PolarPos + PCA_PolarNeg +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_NonPolarPos + PCA_NonPolarNeg + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')


ggsave("Figures/Fig1/PCA_Plots.eps",plot=PCA_plots,device="eps",width=15,height=10,dpi=600)
