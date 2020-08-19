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
colnames(dataNormt) <- rownames(dataNorm.edit) #newly added

#Load vegan
library(vegan)
dist.Metab <- vegdist(dataNormt, method="bray")

#Create groups
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
library("RVAideMemoire")
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.perm.manova(dist.Metab,condI$Cond,nperm=999,p.method="fdr")

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

#####Perform MVS
#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Obtain exometabolite IDs for each strain
BtFinalIDsPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarPos.rds")
CvFinalIDsPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarPos.rds")
PsFinalIDsPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarPos.rds")

#Read in and filter pre-normalized exometabolite table
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarPosMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsPolarPos,CvFinalIDsPolarPos,PsFinalIDsPolarPos))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]

#Match IDs to row number
Btrows <- which(exoMetabolomeReleased$ID %in% BtFinalIDsPolarPos)
Cvrows <- which(exoMetabolomeReleased$ID %in% CvFinalIDsPolarPos)
Psrows <- which(exoMetabolomeReleased$ID %in% PsFinalIDsPolarPos)

#Extract rows and columns of interest for each strain
dataNorm.edit.Bt <- dataNorm.edit[Btrows,1:23]
dataNorm.edit.Cv <- dataNorm.edit[Cvrows,24:47]
dataNorm.edit.Ps <- dataNorm.edit[Psrows,48:71]

#Convert characters to numeric
dataNorm.edit.Bt <- as.data.frame(lapply(dataNorm.edit.Bt, as.numeric))
dataNorm.edit.Cv <- as.data.frame(lapply(dataNorm.edit.Cv, as.numeric))
dataNorm.edit.Ps <- as.data.frame(lapply(dataNorm.edit.Ps, as.numeric))

#Transpose
dataNormt.Bt <- t(dataNorm.edit.Bt)
dataNormt.Cv <- t(dataNorm.edit.Cv)
dataNormt.Ps <- t(dataNorm.edit.Ps)

#Load vegan
library(vegan)
dist.Metab.Bt <- vegdist(dataNormt.Bt, method="bray")
dist.Metab.Cv <- vegdist(dataNormt.Cv, method="bray")
dist.Metab.Ps <- vegdist(dataNormt.Ps, method="bray")

#Extract metaData for each strain
metaD.Bt <- metaD[1:23,]
metaD.Cv <- metaD[24:47,]
metaD.Ps <- metaD[48:71,]

#Calculate betadisper
mod.Bt <- betadisper(dist.Metab.Bt, metaD.Bt)
mod.Cv <- betadisper(dist.Metab.Cv, metaD.Cv)
mod.Ps <- betadisper(dist.Metab.Ps, metaD.Ps)

#Prep for repeated measures anova
metaD.Bt.df <- data.frame(Time=metaD.Bt)
metaD.Cv.df <- data.frame(Time=metaD.Cv)
metaD.Ps.df <- data.frame(Time=metaD.Ps)

#Test for homogeneity of multivariate dispersions
permutest(mod.Bt)
permutest(mod.Cv)
permutest(mod.Ps)

#Null is not rejected for any strain

#Run permutational anova
Bt.permAnova <- adonis2(formula = dist.Metab.Bt ~ Time, data = metaD.Bt.df, permutations = 999, method = "bray")
Cv.permAnova <- adonis2(formula = dist.Metab.Cv ~ Time, data = metaD.Cv.df, permutations = 999, method = "bray")
Ps.permAnova <- adonis2(formula = dist.Metab.Ps ~ Time, data = metaD.Ps.df, permutations = 999, method = "bray")

library("RVAideMemoire")
pairwise.perm.manova(dist.Metab.Bt,metaD.Bt.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Cv,metaD.Cv.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Ps,metaD.Ps.df$Time,nperm=999,p.method="fdr")


#Calculate distances between groups centroids
library(usedist)

#distances between groups centroids all compared to the initial time point
Bt.dist.25toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_25hr.1","X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"))
Bt.dist.30toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"))
Cv.dist.30toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_25hr.1","X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"))
Ps.dist.30toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_45hr.1","X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))

#distances between groups centroids in a timestep manner
Bt.dist.25toInt <- Bt.dist.25toInt
Bt.dist.30to25 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_25hr.1","X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"),
  c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35to30 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"),
  c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40to35 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45to40 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- Cv.dist.25toInt
Cv.dist.30to25 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35to30 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40to35 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45to40 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- Ps.dist.25toInt
Ps.dist.30to25 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_25hr.1","X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"),
  c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35to30 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"),
  c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40to35 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"),
  c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45to40 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"),
  c("X1mem_Ps_45hr.1","X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))


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
library("RVAideMemoire")
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.perm.manova(dist.Metab,condI$Cond,nperm=999,p.method="fdr")

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

#####Perform MVS
#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/PolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Obtain exometabolite IDs for each strain
BtFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
CvFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
PsFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

#Read in and filter pre-normalized exometabolite table
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsPolarNeg,CvFinalIDsPolarNeg,PsFinalIDsPolarNeg))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]

#Match IDs to row number
Btrows <- which(exoMetabolomeReleased$ID %in% BtFinalIDsPolarNeg)
Cvrows <- which(exoMetabolomeReleased$ID %in% CvFinalIDsPolarNeg)
Psrows <- which(exoMetabolomeReleased$ID %in% PsFinalIDsPolarNeg)

#Extract rows and columns of interest for each strain
dataNorm.edit.Bt <- dataNorm.edit[Btrows,1:23]
dataNorm.edit.Cv <- dataNorm.edit[Cvrows,24:47]
dataNorm.edit.Ps <- dataNorm.edit[Psrows,48:71]

#Convert characters to numeric
dataNorm.edit.Bt <- as.data.frame(lapply(dataNorm.edit.Bt, as.numeric))
dataNorm.edit.Cv <- as.data.frame(lapply(dataNorm.edit.Cv, as.numeric))
dataNorm.edit.Ps <- as.data.frame(lapply(dataNorm.edit.Ps, as.numeric))

#Transpose
dataNormt.Bt <- t(dataNorm.edit.Bt)
dataNormt.Cv <- t(dataNorm.edit.Cv)
dataNormt.Ps <- t(dataNorm.edit.Ps)

#Load vegan
library(vegan)
dist.Metab.Bt <- vegdist(dataNormt.Bt, method="bray")
dist.Metab.Cv <- vegdist(dataNormt.Cv, method="bray")
dist.Metab.Ps <- vegdist(dataNormt.Ps, method="bray")

#Extract metaData for each strain
metaD.Bt <- metaD[1:23,]
metaD.Cv <- metaD[24:47,]
metaD.Ps <- metaD[48:71,]

#Calculate betadisper
mod.Bt <- betadisper(dist.Metab.Bt, metaD.Bt)
mod.Cv <- betadisper(dist.Metab.Cv, metaD.Cv)
mod.Ps <- betadisper(dist.Metab.Ps, metaD.Ps)

#Prep for repeated measures anova
metaD.Bt.df <- data.frame(Time=metaD.Bt)
metaD.Cv.df <- data.frame(Time=metaD.Cv)
metaD.Ps.df <- data.frame(Time=metaD.Ps)

#Test for homogeneity of multivariate dispersions
permutest(mod.Bt)
permutest(mod.Cv)
permutest(mod.Ps)

#Null is not rejected for any strain

#Run permutational anova
Bt.permAnova <- adonis2(formula = dist.Metab.Bt ~ Time, data = metaD.Bt.df, permutations = 999, method = "bray")
Cv.permAnova <- adonis2(formula = dist.Metab.Cv ~ Time, data = metaD.Cv.df, permutations = 999, method = "bray")
Ps.permAnova <- adonis2(formula = dist.Metab.Ps ~ Time, data = metaD.Ps.df, permutations = 999, method = "bray")

library("RVAideMemoire")
pairwise.perm.manova(dist.Metab.Bt,metaD.Bt.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Cv,metaD.Cv.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Ps,metaD.Ps.df$Time,nperm=999,p.method="fdr")

#Calculate distances between groups centroids
library(usedist)

#distances between groups centroids all compared to the initial time point
Bt.dist.25toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_25hr.1","X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"))
Bt.dist.30toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.1","X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"))
Cv.dist.30toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.1","X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_25hr.1","X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"))
Ps.dist.30toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.1","X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3","X1mem_Ps_12pt5hr.4"),
  c("X1mem_Ps_45hr.1","X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))

#distances between groups centroids in a timestep manner
Bt.dist.25toInt <- Bt.dist.25toInt
Bt.dist.30to25 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_25hr.1","X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"),
  c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35to30 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_30hr.1","X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"),
  c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40to35 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_35hr.1","X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45to40 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- Cv.dist.25toInt
Cv.dist.30to25 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35to30 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40to35 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45to40 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- Ps.dist.25toInt
Ps.dist.30to25 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_25hr.1","X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"),
  c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35to30 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_30hr.1","X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"),
  c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40to35 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_35hr.1","X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"),
  c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45to40 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_40hr.1","X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"),
  c("X1mem_Ps_45hr.1","X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))




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
library("RVAideMemoire")
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.perm.manova(dist.Metab,condI$Cond,nperm=999,p.method="fdr")

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

#####Perform MVS
#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarPos_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Obtain exometabolite IDs for each strain
BtFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarPos.rds")
CvFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarPos.rds")
PsFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarPos.rds")

#Read in and filter pre-normalized exometabolite table
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarPosMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsNonPolarPos,CvFinalIDsNonPolarPos,PsFinalIDsNonPolarPos))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]

#Match IDs to row number
Btrows <- which(exoMetabolomeReleased$ID %in% BtFinalIDsNonPolarPos)
Cvrows <- which(exoMetabolomeReleased$ID %in% CvFinalIDsNonPolarPos)
Psrows <- which(exoMetabolomeReleased$ID %in% PsFinalIDsNonPolarPos)

#Extract rows and columns of interest for each strain
dataNorm.edit.Bt <- dataNorm.edit[Btrows,1:19]
dataNorm.edit.Cv <- dataNorm.edit[Cvrows,20:42]
dataNorm.edit.Ps <- dataNorm.edit[Psrows,43:59]

#Convert characters to numeric
dataNorm.edit.Bt <- as.data.frame(lapply(dataNorm.edit.Bt, as.numeric))
dataNorm.edit.Cv <- as.data.frame(lapply(dataNorm.edit.Cv, as.numeric))
dataNorm.edit.Ps <- as.data.frame(lapply(dataNorm.edit.Ps, as.numeric))

#Transpose
dataNormt.Bt <- t(dataNorm.edit.Bt)
dataNormt.Cv <- t(dataNorm.edit.Cv)
dataNormt.Ps <- t(dataNorm.edit.Ps)

#Load vegan
library(vegan)
dist.Metab.Bt <- vegdist(dataNormt.Bt, method="bray")
dist.Metab.Cv <- vegdist(dataNormt.Cv, method="bray")
dist.Metab.Ps <- vegdist(dataNormt.Ps, method="bray")

#Extract metaData for each strain
metaD.Bt <- metaD[1:19,]
metaD.Cv <- metaD[20:42,]
metaD.Ps <- metaD[43:59,]

#Calculate betadisper
mod.Bt <- betadisper(dist.Metab.Bt, metaD.Bt)
mod.Cv <- betadisper(dist.Metab.Cv, metaD.Cv)
mod.Ps <- betadisper(dist.Metab.Ps, metaD.Ps)

#Prep for repeated measures anova
metaD.Bt.df <- data.frame(Time=metaD.Bt)
metaD.Cv.df <- data.frame(Time=metaD.Cv)
metaD.Ps.df <- data.frame(Time=metaD.Ps)

#Test for homogeneity of multivariate dispersions
permutest(mod.Bt)
permutest(mod.Cv)
permutest(mod.Ps)

#Null is not rejected for any strain

#Run permutational anova
Bt.permAnova <- adonis2(formula = dist.Metab.Bt ~ Time, data = metaD.Bt.df, permutations = 999, method = "bray")
Cv.permAnova <- adonis2(formula = dist.Metab.Cv ~ Time, data = metaD.Cv.df, permutations = 999, method = "bray")
Ps.permAnova <- adonis2(formula = dist.Metab.Ps ~ Time, data = metaD.Ps.df, permutations = 999, method = "bray")

library("RVAideMemoire")
pairwise.perm.manova(dist.Metab.Bt,metaD.Bt.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Cv,metaD.Cv.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Ps,metaD.Ps.df$Time,nperm=999,p.method="fdr")

#Calculate distances between groups centroids
library(usedist)

#distances between groups centroids all compared to the initial time point
Bt.dist.25toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"))
Bt.dist.30toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"))
Cv.dist.30toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"))
Ps.dist.30toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))

#distances between groups centroids in a timestep manner
Bt.dist.25toInt <- Bt.dist.25toInt
Bt.dist.30to25 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"),
  c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35to30 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"),
  c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40to35 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"),
  c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45to40 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_40hr.1","X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"),
  c("X1mem_Bt_45hr.1","X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- Cv.dist.25toInt
Cv.dist.30to25 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_25hr.1","X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35to30 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40to35 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"),
  c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45to40 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_40hr.1","X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"),
  c("X1mem_Cv_45hr.1","X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- Ps.dist.25toInt
Ps.dist.30to25 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"),
  c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35to30 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"),
  c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40to35 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"),
  c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45to40 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"),
  c("X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))




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
library("RVAideMemoire")
#condI <- unite(condI, newcol, c(Condition, Time), remove=FALSE)
pairwise.perm.manova(dist.Metab,condI$Cond,nperm=999,p.method="fdr")

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

#####Perform MVS
#Read in normalized data from Metaboanalyst
dataNorm <- read.csv("MassSpec/releaseAnalysis/MS/outputFiles/manualEdits/NonPolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)

#Remove first row
dataNorm.edit <- dataNorm[-1,]

#Obtain exometabolite IDs for each strain
BtFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarNeg.rds")
CvFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarNeg.rds")
PsFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarNeg.rds")

#Read in and filter pre-normalized exometabolite table
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsNonPolarNeg,CvFinalIDsNonPolarNeg,PsFinalIDsNonPolarNeg))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]

#Match IDs to row number
Btrows <- which(exoMetabolomeReleased$ID %in% BtFinalIDsNonPolarNeg)
Cvrows <- which(exoMetabolomeReleased$ID %in% CvFinalIDsNonPolarNeg)
Psrows <- which(exoMetabolomeReleased$ID %in% PsFinalIDsNonPolarNeg)

#Extract rows and columns of interest for each strain
dataNorm.edit.Bt <- dataNorm.edit[Btrows,1:17]
dataNorm.edit.Cv <- dataNorm.edit[Cvrows,18:37]
dataNorm.edit.Ps <- dataNorm.edit[Psrows,38:54]

#Convert characters to numeric
dataNorm.edit.Bt <- as.data.frame(lapply(dataNorm.edit.Bt, as.numeric))
dataNorm.edit.Cv <- as.data.frame(lapply(dataNorm.edit.Cv, as.numeric))
dataNorm.edit.Ps <- as.data.frame(lapply(dataNorm.edit.Ps, as.numeric))

#Transpose
dataNormt.Bt <- t(dataNorm.edit.Bt)
dataNormt.Cv <- t(dataNorm.edit.Cv)
dataNormt.Ps <- t(dataNorm.edit.Ps)

#Load vegan
library(vegan)
dist.Metab.Bt <- vegdist(dataNormt.Bt, method="bray")
dist.Metab.Cv <- vegdist(dataNormt.Cv, method="bray")
dist.Metab.Ps <- vegdist(dataNormt.Ps, method="bray")

#Extract metaData for each strain
metaD.Bt <- metaD[1:17,]
metaD.Cv <- metaD[18:37,]
metaD.Ps <- metaD[38:54,]

#Calculate betadisper
mod.Bt <- betadisper(dist.Metab.Bt, metaD.Bt)
mod.Cv <- betadisper(dist.Metab.Cv, metaD.Cv)
mod.Ps <- betadisper(dist.Metab.Ps, metaD.Ps)

#Prep for repeated measures anova
metaD.Bt.df <- data.frame(Time=metaD.Bt)
metaD.Cv.df <- data.frame(Time=metaD.Cv)
metaD.Ps.df <- data.frame(Time=metaD.Ps)

#Test for homogeneity of multivariate dispersions
permutest(mod.Bt)
permutest(mod.Cv)
permutest(mod.Ps)

#Null is not rejected for any strain

#Run permutational anova
Bt.permAnova <- adonis2(formula = dist.Metab.Bt ~ Time, data = metaD.Bt.df, permutations = 999, method = "bray")
Cv.permAnova <- adonis2(formula = dist.Metab.Cv ~ Time, data = metaD.Cv.df, permutations = 999, method = "bray")
Ps.permAnova <- adonis2(formula = dist.Metab.Ps ~ Time, data = metaD.Ps.df, permutations = 999, method = "bray")

library("RVAideMemoire")
pairwise.perm.manova(dist.Metab.Bt,metaD.Bt.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Cv,metaD.Cv.df$Time,nperm=999,p.method="fdr")
pairwise.perm.manova(dist.Metab.Ps,metaD.Ps.df$Time,nperm=999,p.method="fdr")

#Calculate distances between groups centroids
library(usedist)

#distances between groups centroids all compared to the initial time point
Bt.dist.25toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"))
Bt.dist.30toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45toInt <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_12pt5hr.2","X1mem_Bt_12pt5hr.3","X1mem_Bt_12pt5hr.4"),
  c("X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"))
Cv.dist.30toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45toInt <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_12pt5hr.2","X1mem_Cv_12pt5hr.3","X1mem_Cv_12pt5hr.4"),
  c("X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"))
Ps.dist.30toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45toInt <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_12pt5hr.2","X1mem_Ps_12pt5hr.3"),
  c("X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))

#distances between groups centroids in a timestep manner
Bt.dist.25toInt <- Bt.dist.25toInt
Bt.dist.30to25 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_25hr.2","X1mem_Bt_25hr.3","X1mem_Bt_25hr.4"),
  c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"))
Bt.dist.35to30 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_30hr.2","X1mem_Bt_30hr.3","X1mem_Bt_30hr.4"),
  c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"))
Bt.dist.40to35 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_35hr.2","X1mem_Bt_35hr.3","X1mem_Bt_35hr.4"),
  c("X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"))
Bt.dist.45to40 <- dist_between_centroids(dist.Metab.Bt,c("X1mem_Bt_40hr.2","X1mem_Bt_40hr.3","X1mem_Bt_40hr.4"),
  c("X1mem_Bt_45hr.2","X1mem_Bt_45hr.3"))

Cv.dist.25toInt <- Cv.dist.25toInt
Cv.dist.30to25 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_25hr.2","X1mem_Cv_25hr.3","X1mem_Cv_25hr.4"),
  c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"))
Cv.dist.35to30 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_30hr.1","X1mem_Cv_30hr.2","X1mem_Cv_30hr.3","X1mem_Cv_30hr.4"),
  c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"))
Cv.dist.40to35 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_35hr.1","X1mem_Cv_35hr.2","X1mem_Cv_35hr.3","X1mem_Cv_35hr.4"),
  c("X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"))
Cv.dist.45to40 <- dist_between_centroids(dist.Metab.Cv,c("X1mem_Cv_40hr.2","X1mem_Cv_40hr.3","X1mem_Cv_40hr.4"),
  c("X1mem_Cv_45hr.2","X1mem_Cv_45hr.3","X1mem_Cv_45hr.4"))

Ps.dist.25toInt <- Ps.dist.25toInt
Ps.dist.30to25 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_25hr.2","X1mem_Ps_25hr.3","X1mem_Ps_25hr.4"),
  c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"))
Ps.dist.35to30 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_30hr.2","X1mem_Ps_30hr.3","X1mem_Ps_30hr.4"),
  c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"))
Ps.dist.40to35 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_35hr.2","X1mem_Ps_35hr.3","X1mem_Ps_35hr.4"),
  c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"))
Ps.dist.45to40 <- dist_between_centroids(dist.Metab.Ps,c("X1mem_Ps_40hr.2","X1mem_Ps_40hr.3","X1mem_Ps_40hr.4"),
  c("X1mem_Ps_45hr.2","X1mem_Ps_45hr.3","X1mem_Ps_45hr.4"))




#Let's put all mass spec analyses together

library(patchwork)

PCA_plots <- PCA_PolarPos + PCA_PolarNeg +
  theme(legend.position="right",legend.text = element_text(size = 20)) + guides(color=guide_legend(title = "Strain",order=1),fill=FALSE,size=guide_legend(title = "Time(h)",order=3)) +
  PCA_NonPolarPos + PCA_NonPolarNeg + plot_layout(ncol=2,nrow=2) + plot_annotation(tag_levels = 'A')


ggsave("Figures/Fig1/PCA_Plots.eps",plot=PCA_plots,device="eps",width=15,height=10,dpi=600)
