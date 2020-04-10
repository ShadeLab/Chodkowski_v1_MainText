#####This is analysis is determine the number of released exometabolites in each isolate#####
#####and, determine the identity of a subset of these exometabolites#####
#####This is for Nonpolar analysis in negative ionization mode#####

library(dplyr)

#Read in MZmine final feature table for MS analysis
#####################MAKE SURE YOU CHANGE THIS ON GITHUB FOR FINAL CODE SUBMISSION#####################
filt <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegwExConForQC.csv",check.names=F)

#Remove outlier-1mem_Bt_45hr-4
filt_f <- filt[,-c(119)]

#Retain only monoculture samples
library(reshape)
filtMelt <- melt(filt_f,id.vars=c("ID","RT","MZ"))

#Split column name string so we get 1mem, 2mem, 3mem and ExCon identifiers
msD <- data.frame(do.call('rbind', strsplit(as.character(filtMelt$variable),'_',fixed=TRUE)))
#Place this into melted object
filtMelt$group <- msD$X1

#Retain only 1mem and ExControl in melted object
keep <- c("1mem","ExControl-ExControl")
filtMeltMono <- filtMelt[which(filtMelt$group %in% keep),]

#Remove group identifier
remove <- "group"
filtMeltMonoFinal <- filtMeltMono[,-which(names(filtMeltMono) %in% remove)]

#Convert back to dataFrame and prepare monoculture files
library(reshape2)
exoMetabMonoCulture <- dcast(filtMeltMonoFinal, ID + MZ + RT ~ variable)
write.csv(exoMetabMonoCulture,"MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",row.names=FALSE)
exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",check.names=F)
#Order columns
exoMetabolomeMonoCultureOrdered <- exoMetabolomeMonoCulture[ , order(names(exoMetabolomeMonoCulture))]
write.csv(exoMetabolomeMonoCultureOrdered,"MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonocultureOrdered.csv",row.names=FALSE)

#####Start of released exometabolite analysis#####
exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",check.names=F)

library(reshape)
#Melt dataframe
filtM <- melt(exoMetabMonoCulture,id.vars=c("ID","RT","MZ"))

#Change NA values and 0s to low value
min(filtM[filtM$value>0,5],na.rm=TRUE)
filtM[,5][filtM[,5]==0] <-1/1000*min(filtM[filtM$value>0,5],na.rm=TRUE)

###################################Filter 1###################################

#Remove features where the maximum value is found in the external control

#Create groups for averaging
msD <- data.frame(do.call('rbind', strsplit(as.character(filtM$variable),'-',fixed=TRUE)))
#Place this into melted object
filtM$groupTime <- msD$X1

#Obtain the group with the maximum peak area value
FeatureMax <- filtM %>%
             group_by(ID) %>%
             filter(value == max(value,na.rm=T))

#Obtain feature IDs where the max value was found in the external control

MaxInCon <- FeatureMax[FeatureMax$groupTime=="ExControl",]

#Remove these features
filtS1 <- exoMetabMonoCulture[-c(which(exoMetabMonoCulture$ID %in% MaxInCon$ID)),]
#6411 features remain

###################################Filter 2###################################

#Keep features that have a positive correlation with time

#Melt dataframe
filtMelt <- melt(filtS1,id.vars=c("ID","RT","MZ"))

#Change NA values and 0s to 2- will perform log2 variance stabilization so 2
filtMelt[,5][filtMelt[,5]<2] <- 2

#We need to average first before we log transform
msD <- data.frame(do.call('rbind', strsplit(as.character(filtMelt$variable),'-',fixed=TRUE)))
#Place this into melted object
filtMelt$groupTime <- msD$X1

#Remove ExControls
filtMeltCon<- filtMelt[-c(which(filtMelt$groupTime=="ExControl")),]
#Remove varible column
filtMeltCon$variable <- NULL

library(dplyr)
detach("package:plyr", unload=TRUE)
groupedMono <- filtMeltCon %>% group_by(ID,RT,MZ,groupTime)
#Get mean and std
groupedMonoSum <- summarise(groupedMono, mean=mean(value))

#Add conditions
#each 6 refers to 6 time points while the second number is the total features remaining from filtering S1
groupedMonoSum$conditions <- rep(c("Bt","Cv","Ps"),each=6,6411)
#Add time
groupedMonoSum$time <- rep(c(12.5,25,30,35,40,45),19233) #This number is the number above * 3 (# of conditions)

#variance stablization
groupedMonoSum$mean <- log2(groupedMonoSum$mean)

#Extract individual monocultures
groupedMonoSumBt <- groupedMonoSum[groupedMonoSum$conditions=="Bt" ,]
groupedMonoSumCv<- groupedMonoSum[groupedMonoSum$conditions=="Cv",]
groupedMonoSumPs<- groupedMonoSum[groupedMonoSum$conditions=="Ps",]

#Test for correlations
groupedMonoBt <- groupedMonoSumBt %>% group_by(ID)
groupedMonoCv <- groupedMonoSumCv %>% group_by(ID)
groupedMonoPs <- groupedMonoSumPs %>% group_by(ID)

#Get mean and std
#Detach plyr if not working
groupedMonoBtcorr <- summarize(groupedMonoBt, cor(mean,time))
groupedMonoCvcorr <- summarize(groupedMonoCv, cor(mean,time))
groupedMonoPscorr <- summarize(groupedMonoPs, cor(mean,time))

#Get positive correlated exometabolites
BtIDS <- groupedMonoBtcorr[groupedMonoBtcorr$"cor(mean, time)">0,]$ID
CvIDS <- groupedMonoCvcorr[groupedMonoCvcorr$"cor(mean, time)">0,]$ID
PsIDS <- groupedMonoPscorr[groupedMonoPscorr$"cor(mean, time)">0,]$ID
#Contains NAs because the mean was the same value across all time points. E.g.
#cor(c(1,1,1,1), c(1,2,3,4))
#[1] NA
#Warning message:
#In cor(c(1, 1, 1, 1), c(1, 2, 3, 4)) : the standard deviation is zero
#We can ignore this warning.

#Obtain IDs with pos cor exometabolites
groupedMonoSumBtCorrF <- groupedMonoSumBt[c(which(groupedMonoSumBt$ID %in% BtIDS)),]
groupedMonoSumCvCorrF <- groupedMonoSumCv[c(which(groupedMonoSumCv$ID %in% CvIDS)),]
groupedMonoSumPsCorrF <- groupedMonoSumPs[c(which(groupedMonoSumPs$ID %in% PsIDS)),]

#Make sure that you have the correct number of filtered rows. groupedMonoBtcorr (or Cv/Ps) have a lot of NAs with this object
#The true number of features left should be:
#length(unique(BtIDS)) - 1 = nrow(groupedMonoSumBtCorrF)/6 #Looks good
#length(unique(CvIDS)) - 1 = nrow(groupedMonoSumCvCorrF)/6 #Looks good
#length(unique(PsIDS)) - 1 = nrow(groupedMonoSumPsCorrF)/6 #Looks good

###################################Filter 3###################################

#Make sure the first time point is the minimum value. We're looking for compounds released during stationary phase

FeatureMinBt <- groupedMonoSumBtCorrF %>%
             group_by(ID) %>%
             filter(mean == min(mean,na.rm=T))

FeatureMinBtExpo <- FeatureMinBt[FeatureMinBt$time=="12.5",]

FeatureMinCv <- groupedMonoSumCvCorrF %>%
             group_by(ID) %>%
             filter(mean == min(mean,na.rm=T))

FeatureMinCvExpo <- FeatureMinCv[FeatureMinCv$time=="12.5",]

FeatureMinPs <- groupedMonoSumPsCorrF %>%
             group_by(ID) %>%
             filter(mean == min(mean,na.rm=T))

FeatureMinPsExpo <- FeatureMinPs[FeatureMinPs$time=="12.5",]

groupedMonoBtminExpo <- groupedMonoSumBtCorrF[c(which(groupedMonoSumBtCorrF$ID %in% FeatureMinBtExpo$ID)),]
groupedMonoCvminExpo <- groupedMonoSumCvCorrF[c(which(groupedMonoSumCvCorrF$ID %in% FeatureMinCvExpo$ID)),]
groupedMonoPsminExpo <- groupedMonoSumPsCorrF[c(which(groupedMonoSumPsCorrF$ID %in% FeatureMinPsExpo$ID)),]

###################################Filter 4###################################

#LogFoldChange (LFC) filter. LFC from last timepoint - first timepoint should be > 1.

BtmonoLFC <- groupedMonoBtminExpo[groupedMonoBtminExpo$time==c("12.5","45"),]
CvmonoLFC <- groupedMonoCvminExpo[groupedMonoCvminExpo$time==c("12.5","45"),]
PsmonoLFC <- groupedMonoPsminExpo[groupedMonoPsminExpo$time==c("12.5","45"),]

BtmonoLFCDiff <- BtmonoLFC %>%
             group_by(ID) %>%
             mutate(
             Difference = mean - lag(mean)
             )
table(BtmonoLFCDiff$Difference>1)

CvmonoLFCDiff <- CvmonoLFC %>%
             group_by(ID) %>%
             mutate(
             Difference = mean - lag(mean)
             )
table(CvmonoLFCDiff$Difference>1)

PsmonoLFCDiff <- PsmonoLFC %>%
             group_by(ID) %>%
             mutate(
             Difference = mean - lag(mean)
             )
table(PsmonoLFCDiff$Difference>1)

#Obtain features that have a LFC > 1
BtFeat <- groupedMonoBtminExpo[which(groupedMonoBtminExpo$ID %in% unique(BtmonoLFCDiff$ID[BtmonoLFCDiff$Difference>1])),]
CvFeat <- groupedMonoCvminExpo[which(groupedMonoCvminExpo$ID %in% unique(CvmonoLFCDiff$ID[CvmonoLFCDiff$Difference>1])),]
PsFeat <- groupedMonoPsminExpo[which(groupedMonoPsminExpo$ID %in% unique(PsmonoLFCDiff$ID[PsmonoLFCDiff$Difference>1])),]

###################################Filter 5###################################

#Exometabolites in stationary phase have correlations >0.7

#Remove the exponential phase time point
BtFeatSP <- BtFeat[which(BtFeat$time!="12.5"),]
CvFeatSP <- CvFeat[which(CvFeat$time!="12.5"),]
PsFeatSP <- PsFeat[which(PsFeat$time!="12.5"),]

#Test for correlations
BtFeatSPCor <- BtFeatSP %>% group_by(ID)
CvFeatSPCor <- CvFeatSP %>% group_by(ID)
PsFeatSPCor <- PsFeatSP %>% group_by(ID)

#Get mean and std
#Detach plyr if not working
detach("package:plyr", unload=TRUE)
BtFeatSPCorR <- summarize(BtFeatSPCor, cor(mean,time))
CvFeatSPCorR <- summarize(CvFeatSPCor, cor(mean,time))
PsFeatSPCorR <- summarize(PsFeatSPCor, cor(mean,time))

#Get positive correlated exometabolites
BtcorrFeat <- BtFeatSPCorR[BtFeatSPCorR$"cor(mean, time)">0.7,]$ID
CvcorrFeat <- CvFeatSPCorR[CvFeatSPCorR$"cor(mean, time)">0.7,]$ID
PscorrFeat <- PsFeatSPCorR[PsFeatSPCorR$"cor(mean, time)">0.7,]$ID

#Filter features with correlations greater than 0.7
BtFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% BtcorrFeat),]
CvFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% CvcorrFeat),]
PsFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% PscorrFeat),]

###################################Filter 6###################################

#Noise filter: The min value at TP 45 (last timepoint) should be 3x the max value found in the ExCon (Filter suggested by JGI)

library(reshape)
filtBt <- melt(BtFeatures,id.vars=c("ID","RT","MZ"))
filtCv <- melt(CvFeatures,id.vars=c("ID","RT","MZ"))
filtPs <- melt(PsFeatures,id.vars=c("ID","RT","MZ"))

#Change NA values and 0s to low value
min(filtBt[filtBt$value>0,5],na.rm=TRUE)
min(filtCv[filtCv$value>0,5],na.rm=TRUE)
min(filtPs[filtPs$value>0,5],na.rm=TRUE)

#filtM[,5][is.zero(filtM[,5])] <-1/1000*min(filtM[filtM$value>0,5],na.rm=TRUE)
filtBt[,5][filtBt[,5]==0] <-1/1000*min(filtBt[filtBt$value>0,5],na.rm=TRUE)
filtCv[,5][filtCv[,5]==0] <-1/1000*min(filtCv[filtCv$value>0,5],na.rm=TRUE)
filtPs[,5][filtPs[,5]==0] <-1/1000*min(filtPs[filtPs$value>0,5],na.rm=TRUE)

#Create groups for averaging
msBt <- data.frame(do.call('rbind', strsplit(as.character(filtBt$variable),'-',fixed=TRUE)))
msCv <- data.frame(do.call('rbind', strsplit(as.character(filtCv$variable),'-',fixed=TRUE)))
msPs <- data.frame(do.call('rbind', strsplit(as.character(filtPs$variable),'-',fixed=TRUE)))

#Place this into melted object
filtBt$groupTime <- msBt$X1
filtCv$groupTime <- msCv$X1
filtPs$groupTime <- msPs$X1

#Filter out only 45hr TP for monoculture and ExCon
#Detach plyr if not working
detach("package:plyr", unload=TRUE)
require(dplyr)
Btsamps <- c("1mem_Bt_45hr", "ExControl")
Cvsamps <- c("1mem_Cv_45hr", "ExControl")
Pssamps <- c("1mem_Ps_45hr", "ExControl")

filtMonoConBt <- filter(filtBt, groupTime %in% Btsamps)
filtMonoConCv <- filter(filtCv, groupTime %in% Cvsamps)
filtMonoConPs <- filter(filtPs, groupTime %in% Pssamps)

#Obtain max values for each ID for
require(dplyr)
filtMonoConBtMax <- filtMonoConBt %>% group_by(ID,groupTime) %>% summarise(value = max(value))
filtMonoConCvMax <- filtMonoConCv %>% group_by(ID,groupTime) %>% summarise(value = max(value))
filtMonoConPsMax <- filtMonoConPs %>% group_by(ID,groupTime) %>% summarise(value = max(value))

#Obtain min values for each ID for
filtMonoConBtMin <- filtMonoConBt %>% group_by(ID,groupTime) %>% summarise(value = min(value))
filtMonoConCvMin <- filtMonoConCv %>% group_by(ID,groupTime) %>% summarise(value = min(value))
filtMonoConPsMin <- filtMonoConPs %>% group_by(ID,groupTime) %>% summarise(value = min(value))

#Obtain and combine min from 45hr and max from excontrol
filtMonoConBtMaxMin <- rbind(filtMonoConBtMin[filtMonoConBtMin$groupTime=="1mem_Bt_45hr",],filtMonoConBtMax[filtMonoConBtMax$groupTime=="ExControl",])
filtMonoConCvMaxMin <- rbind(filtMonoConCvMin[filtMonoConCvMin$groupTime=="1mem_Cv_45hr",],filtMonoConCvMax[filtMonoConCvMax$groupTime=="ExControl",])
filtMonoConPsMaxMin <- rbind(filtMonoConPsMin[filtMonoConPsMin$groupTime=="1mem_Ps_45hr",],filtMonoConPsMax[filtMonoConPsMax$groupTime=="ExControl",])

#Obtain fold change of max values between 45 TP and ExCon
BtmonoFC <- filtMonoConBtMaxMin %>%
             group_by(ID) %>%
             mutate(
             FC = lag(value)/value
             )
table(BtmonoFC$FC>3)

CvmonoFC <- filtMonoConCvMaxMin %>%
             group_by(ID) %>%
             mutate(
             FC = lag(value)/value
             )
table(CvmonoFC$FC>3)

PsmonoFC <- filtMonoConPsMaxMin %>%
             group_by(ID) %>%
             mutate(
             FC = lag(value)/value
             )
table(PsmonoFC$FC>3)

#Doublecheck all is well
BtmonoFC[order(BtmonoFC$ID),]
CvmonoFC[order(CvmonoFC$ID),]
PsmonoFC[order(PsmonoFC$ID),]

BtMinMaxIDs <- unique(BtmonoFC$ID[BtmonoFC$FC>3]) #1585
CvMinMaxIDs <- unique(CvmonoFC$ID[CvmonoFC$FC>3]) #358
PsMinMaxIDs <- unique(PsmonoFC$ID[PsmonoFC$FC>3]) #348

BtMinMaxFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% BtMinMaxIDs),]
CvMinMaxFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% CvMinMaxIDs),]
PsMinMaxFeatures <- exoMetabMonoCulture[which(exoMetabMonoCulture$ID %in% PsMinMaxIDs),]

###################################Filter 7###################################

#RSD filter below 20%- All time points needs to be below 20%

filtMeltBT <- melt(BtMinMaxFeatures,id.vars=c("ID","RT","MZ"))
filtMeltCV <- melt(CvMinMaxFeatures,id.vars=c("ID","RT","MZ"))
filtMeltPS <- melt(PsMinMaxFeatures,id.vars=c("ID","RT","MZ"))

#Change NA values and 0s to low value
filtMeltBT[,5][filtMeltBT[,5]<2] <- 2
filtMeltCV[,5][filtMeltCV[,5]<2] <- 2
filtMeltPS[,5][filtMeltPS[,5]<2] <- 2

#Perform RSD filter
rsdBT <- filtMeltBT
rsdCV <- filtMeltCV
rsdPS <- filtMeltPS

#log transform values
rsdBT$value <- log2(rsdBT$value)
rsdCV$value <- log2(rsdCV$value)
rsdPS$value <- log2(rsdPS$value)

#Create tibble
msDBT <- data.frame(do.call('rbind', strsplit(as.character(rsdBT$variable),'-',fixed=TRUE)))
msDCV <- data.frame(do.call('rbind', strsplit(as.character(rsdCV$variable),'-',fixed=TRUE)))
msDPS <- data.frame(do.call('rbind', strsplit(as.character(rsdPS$variable),'-',fixed=TRUE)))

#Place this into melted object
rsdBT$groupTime <- msDBT$X1
rsdCV$groupTime <- msDCV$X1
rsdPS$groupTime <- msDPS$X1

#Add conditions to tibble
condBT <- data.frame(do.call('rbind', strsplit(as.character(rsdBT$variable),'_',fixed=TRUE)))
condCV <- data.frame(do.call('rbind', strsplit(as.character(rsdCV$variable),'_',fixed=TRUE)))
condPS <- data.frame(do.call('rbind', strsplit(as.character(rsdPS$variable),'_',fixed=TRUE)))

rsdBT$cond <- condBT$X2
rsdCV$cond <- condCV$X2
rsdPS$cond <- condPS$X2

#Keep only monocultures
rsdBtOnly <- rsdBT[rsdBT$cond=="Bt",]
rsdCvOnly <- rsdCV[rsdCV$cond=="Cv",]
rsdPsOnly <- rsdPS[rsdPS$cond=="Ps",]

library(dplyr)
detach("package:plyr", unload=TRUE)
groupedRSDBt <- rsdBtOnly %>% group_by(ID,RT,MZ,groupTime)
groupedRSDCv <- rsdCvOnly %>% group_by(ID,RT,MZ,groupTime)
groupedRSDPs <- rsdPsOnly %>% group_by(ID,RT,MZ,groupTime)

#Get mean and std
groupRSDSumBt <- summarise(groupedRSDBt, mean=mean(value), sd=sd(value))
groupRSDSumCv <- summarise(groupedRSDCv, mean=mean(value), sd=sd(value))
groupRSDSumPs <- summarise(groupedRSDPs, mean=mean(value), sd=sd(value))

#Calculate RSD
groupRSDSumBt$RSD <-groupRSDSumBt$sd/groupRSDSumBt$mean
groupRSDSumCv$RSD <-groupRSDSumCv$sd/groupRSDSumCv$mean
groupRSDSumPs$RSD <-groupRSDSumPs$sd/groupRSDSumPs$mean

#Determine which IDs have a RSD of 0.
groupRSDSumBT0s <- as.data.frame(groupRSDSumBt[groupRSDSumBt$RSD==0,])
groupRSDSumCv0s <- as.data.frame(groupRSDSumCv[groupRSDSumCv$RSD==0,])
groupRSDSumPs0s <- as.data.frame(groupRSDSumPs[groupRSDSumPs$RSD==0,])

#The plan: get rid of highly variable features. Features across all time point should have a RSD <20%. Standard in LCMS is 20%
#Find max RSD value
groupRSDSumBtMAX <- groupRSDSumBt %>%
             group_by(ID) %>%
             filter(RSD == max(RSD,na.rm=T))

groupRSDSumCvMAX <- groupRSDSumCv %>%
             group_by(ID) %>%
             filter(RSD == max(RSD,na.rm=T))

groupRSDSumPsMAX <- groupRSDSumPs %>%
             group_by(ID) %>%
             filter(RSD == max(RSD,na.rm=T))

#Remove features whose minimum RSD value is >0.2. High variable features
BtFinalIDs <- groupRSDSumBtMAX$ID[groupRSDSumBtMAX$RSD<0.2] #1066
CvFinalIDs <- groupRSDSumCvMAX$ID[groupRSDSumCvMAX$RSD<0.2] #202
PsFinalIDs <- groupRSDSumPsMAX$ID[groupRSDSumPsMAX$RSD<0.2] #221

#Confirm filtered IDs that also had a RSD value of 0 was the RSD value at the initial timepoint
groupRSDSumBT0s[which(groupRSDSumBT0s$ID %in% BtFinalIDs),] #Remove IDs 8516,34428,34545,39976,51530
groupRSDSumCv0s[which(groupRSDSumCv0s$ID %in% CvFinalIDs),] #Looks good
groupRSDSumPs0s[which(groupRSDSumPs0s$ID %in% PsFinalIDs),] #Looks good

BtFinalIDs <- BtFinalIDs[!BtFinalIDs %in% c(8516,34428,34545,39976,51530)] #1061

#Save IDs for each isolate that have been determined to be released
saveRDS(BtFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarNeg.rds")
saveRDS(CvFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarNeg.rds")
saveRDS(PsFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarNeg.rds")

##########End of filtering steps to determine released exometabolites##########





###################Start here to continue with identification###################
BtFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarNeg.rds")
CvFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarNeg.rds")
PsFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarNeg.rds")

#Create a venn diagram
finalFeats <- list(Bt=BtFinalIDsNonPolarNeg,Cv=CvFinalIDsNonPolarNeg,Ps=PsFinalIDsNonPolarNeg)
finalFeatsV <- lapply(finalFeats, function(x) x[!is.na(x)])
library(VennDiagram)

#Supplementary Figure 1
venn.diagram(finalFeatsV,filename="Figures/SFig1/MonoOverlappingFeaturesNonPolarNeg.tiff")

#Venn output shows
#995 unique features for Bt
#151 unique features for Cv
#176 unique features for Ps

#####Use venn diagrams from Polar/NonPolar and Positive/Negative modes (Supplemental Figure 1) to generate Table 1#####

#Let's export a CSV file for all released exometabolites for Metaboanalyst

#Read in file
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsNonPolarNeg,CvFinalIDsNonPolarNeg,PsFinalIDsNonPolarNeg))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]
#Add internal standard- 2­Amino­3­bromo­5­methylbenzoic acid, ID 14
exoMetabolomeReleased <- rbind(exoMetabolomeReleased,exoMetabolome[exoMetabolome$ID=="14",])

#Re-order columns
exoMetabolomeReleasedO <- exoMetabolomeReleased[ , order(names(exoMetabolomeReleased))]

#Make IDs rownames
rownames(exoMetabolomeReleasedO) <- exoMetabolomeReleasedO$ID

#Extract strain to separate dataframes
Bt.ext <- exoMetabolomeReleasedO[, grep('Bt', names(exoMetabolomeReleasedO))]
Cv.ext <- exoMetabolomeReleasedO[, grep('Cv', names(exoMetabolomeReleasedO))]
Ps.ext <- exoMetabolomeReleasedO[, grep('Ps', names(exoMetabolomeReleasedO))]

#Add standard to final IDs
BtFinalIDsNonPolarNeg <- append(BtFinalIDsNonPolarNeg,14)
CvFinalIDsNonPolarNeg <- append(CvFinalIDsNonPolarNeg,14)
PsFinalIDsNonPolarNeg <- append(PsFinalIDsNonPolarNeg,14)

#Obtain exometabolite IDs that did that fit criteria of released for a particular strain
Bt.notR <- which(!rownames(Bt.ext) %in% BtFinalIDsNonPolarNeg)
Cv.notR <- which(!rownames(Cv.ext) %in% CvFinalIDsNonPolarNeg)
Ps.notR <- which(!rownames(Ps.ext) %in% PsFinalIDsNonPolarNeg)

#Convert values that do not match release criteria to NAs
Bt.ext[Bt.notR,] = NA
Cv.ext[Cv.notR,] = NA
Ps.ext[Ps.notR,] = NA

#Put the dataframe back together
msInfo <- c("ID","MZ","RT")
msInfo.cols <- which(colnames(exoMetabolomeReleasedO) %in% msInfo)
exoMetabolomeReleasedO.f <- cbind(Bt.ext,Cv.ext,Ps.ext,exoMetabolomeReleasedO[,msInfo.cols])

NA_Index  <- as.data.frame(which(is.na(exoMetabolomeReleasedO.f), arr.ind=TRUE))

#Write csv and save indicies of NAs
write.csv(exoMetabolomeReleasedO,"MassSpec/releaseAnalysis/MS/outputFiles/MZmineNonPolarNegallReleasedMetabolitesIndividualSamples.csv")
saveRDS(NA_Index,"MassSpec/releaseAnalysis/MS/outputFiles/NonPolarNeg_NA_Index.rds")

#Note: MZmineNonPolarNegallReleasedMetabolitesIndividualSamples.csv file was then manually edited to prepare for Metaboanalyst.
#This was done is two ways:
#1) MZmineNonPolarNegallReleasedMetabolitesIndividualSamples_readyForMetaboAnalyst_wTime.csv for HeatMap analysis.
    #Perform statistical analysis on MetaboAnalyst and generate heatmap averagign bioreps: https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml
    #Note, pseudosamples were added to the data where only two bioreps were present. This was performed because MetaboAnalyst requries at least 3 samples per group.
    #The pseudosample was simply an average of the two bioreps, so it has no effect on the overall results after averaging bioreps
#2) MZmineNonPolarNegallReleasedMetabolitesIndividualSamples_readyForMetaboAnalyst.csv for PCA analysis
    #Data was normalized in metaboanalyst and normalized data was exported to NonPolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv
      #This was then uploaded in R to perform PCA in vegan, generating Figure 1.

#MZmineNonPolarNegallReleasedMetabolitesIndividualSamples.csv manually edited file was uploade to MetaboAnalyst
#This file was used to generate Figure 1 and Supplemental Figure 2
#Follow the Rscript, "PCAs.R", to finish generating Figure 1.
#Supp. Fig. 2 was generated directly in MetaboAnalyst





#####################Identification of released metabolites#####################
BtFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarNeg.rds")
CvFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarNeg.rds")
PsFinalIDsNonPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarNeg.rds")

exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZmineNonPolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsNonPolarNeg,CvFinalIDsNonPolarNeg,PsFinalIDsNonPolarNeg))

filtJGI <- exoMetabolomeMonoCulture[c(which(exoMetabolomeMonoCulture$ID %in% unique(finalFeatsJGI))),]

###Steps for exometabolite identification
# 1) Nist17
# 2) MoNa
#For any overlaps, #1 overrides. For overlaps between 2 & 3, highest cosine simialrity

#Match to nist17 and MONA
#Load MSMS analysis from MZmine
MSMSdat <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZmine_NonPolarNegMSMS_peakArea.csv",check.names=F)

#They match! Filtered features with MSMS data
BtMSMS <- BtFinalIDsNonPolarNeg[which(BtFinalIDsNonPolarNeg %in% MSMSdat$"row ID")] #236/1061 filtered features
CvMSMS <- CvFinalIDsNonPolarNeg[which(CvFinalIDsNonPolarNeg %in% MSMSdat$"row ID")] #70/202 filtered features
PsMSMS <- PsFinalIDsNonPolarNeg[which(PsFinalIDsNonPolarNeg %in% MSMSdat$"row ID")] #43/221 filtered features

#See which IDs match with MSMS data
MSMSmatchesNist17 <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZmineNonPolarNegMSMSnist17Matches.csv",check.names=F)
MSMSmatchesMONA <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/NonPolarNegMSMS_MONAspectraMatches.csv",check.names=F)

#MSMS matches from nist17
BtMSSpecsNist17 <- BtMSMS[which(BtMSMS %in% MSMSmatchesNist17$"cluster index")]
CvMSSpecsNist17 <- CvMSMS[which(CvMSMS %in% MSMSmatchesNist17$"cluster index")]
PsMSSpecsNist17 <- PsMSMS[which(PsMSMS %in% MSMSmatchesNist17$"cluster index")]

MSSpecsNist17F <- unique(c(BtMSSpecsNist17, CvMSSpecsNist17,PsMSSpecsNist17))

#MSMS matches from MONA
BtMSSpecsMONA <- BtMSMS[which(BtMSMS %in% MSMSmatchesMONA$rowID)]
CvMSSpecsMONA <- CvMSMS[which(CvMSMS %in% MSMSmatchesMONA$rowID)]
PsMSSpecsMONA <- PsMSMS[which(PsMSMS %in% MSMSmatchesMONA$rowID)]

MSSpecsMONAF <- unique(c(BtMSSpecsMONA, CvMSSpecsMONA,PsMSSpecsMONA))

#Of those identified by nist17 or MONA, which ones overlap?
MvN <- MSSpecsNist17F[which(MSSpecsNist17F %in% MSSpecsMONAF)]

N17 <- MSMSmatchesNist17[which(MSMSmatchesNist17$"cluster index" %in% MvN),]
MONA <- MSMSmatchesMONA[which(MSMSmatchesMONA$rowID %in% MvN),]

#Order rows
N17o <- N17[order(N17[,10]),]
MONAo <- MONA[with(MONA,order(rowID)),]

#Which identifications to use from Nist17?
monaRemove <- N17o[N17o$MQScore>MONAo$cos,10]
#Which identifications to use from MONA?
nist17Remove <- MONAo[MONAo$cos>N17o$MQScore,1]

#Remove NIST17 and MONA overlaps
MSSpecsMONAFf <- MSSpecsMONAF[!MSSpecsMONAF %in% monaRemove]
MSSpecsNist17Ff <- MSSpecsNist17F[!MSSpecsNist17F %in% nist17Remove]
which(MSSpecsNist17Ff %in% MSSpecsMONAFf) #none. Good.

MONAhits <- MSMSmatchesMONA[which(MSMSmatchesMONA$rowID %in% MSSpecsMONAFf),]
library(stringr)
MONAhits$LibraryID <- str_split_fixed(MONAhits$molecule, "as", 2)[,1]
nist17hits <- MSMSmatchesNist17[which(MSMSmatchesNist17$"cluster index" %in% MSSpecsNist17Ff),]

#Now we need to rid of redundancies within a method
#For some reason, the mzMine procedure picks redundant peaks. We want to keep the ID with the higher cosine score

#Check for redundant compounds

#Nist17 redundants
length(nist17hits$"cluster index")   #5
length(unique(nist17hits$"cluster index"))  #5
length(nist17hits$LibraryID)   #5
length(unique(nist17hits$LibraryID))  #5
#obtain redundant features
dupeNist17 <- nist17hits[,"LibraryID"]
nist17hits[duplicated(dupeNist17) | duplicated(dupeNist17, fromLast=TRUE),]

#MONA redundants
length(MONAhits$rowID) #10
length(unique(MONAhits$rowID)) #10
length(MONAhits$LibraryID) #10
length(unique(MONAhits$LibraryID)) #9
#obtain redundant features
dupeMONA <- MONAhits[,"LibraryID"]
MONAhits[duplicated(dupeMONA) | duplicated(dupeMONA, fromLast=TRUE),]

#Subset columns
MONAhitsSub <- MONAhits[,c(1,6)]
nist17hitsSub <- nist17hits[,c(10,2)]

#Add additional column for source of identification
MONAhitsSub$source <- rep("mona",nrow(MONAhitsSub))
nist17hitsSub$source <- rep("nist17",nrow(nist17hitsSub))

#Rename columns
names(MONAhitsSub) <- c("ID","Label","source")
names(nist17hitsSub) <- c("ID","Label","source")

#Bring it all together
identifiedFinal <- rbind(MONAhitsSub,nist17hitsSub)

#Obtain IDS from all methods
filtFinalMSMS <- exoMetabolomeMonoCulture[c(which(exoMetabolomeMonoCulture$ID %in% unique(identifiedFinal$ID))),]
#Add internal standard- 2­Amino­3­bromo­5­methylbenzoic acid, ID 14
filtFinalMSMS <- rbind(filtFinalMSMS,exoMetabolomeMonoCulture[exoMetabolomeMonoCulture$ID=="14",])
#Re-order columns
filtFinalMSMSO <- filtFinalMSMS[ , order(names(filtFinalMSMS))]
#Write csv for individual samples
write.csv(filtFinalMSMSO,"MassSpec/releaseAnalysis/MSMS/outputFiles/MZmineNonPolarNegReleasedMetaboliteswIDS_IndividualSamples.csv")

#Write compounds to CSV
write.csv(identifiedFinal,"MassSpec/releaseAnalysis/MSMS/outputFiles/NonPolarNegMSMScompounds.csv")

######################################################Manually edit file######################################################

#First, manually remove redundancies within the same identification #Keep the exometabolite with the higher cosine score
MONAhits[duplicated(dupeMONA) | duplicated(dupeMONA, fromLast=TRUE),]
nist17hits[duplicated(dupeNist17) | duplicated(dupeNist17, fromLast=TRUE),]

#Lastly, look through molecules. Is it logical that it should be there? Specifically for MONA and NIST17 hits.
#For example, remove 7-(2-hydroxypropan-2-yl)-1,4a-dimethyl-2,3,4,9,10,10a-hexahydrophenanthrene-1-carboxylic acid, KOBUSONE,
  #(1S,8R,9R)-8-hydroxy-4-(propan-2-ylidene)-10-oxatricyclo[7.2.1.01,5]dodecane-8-carboxylic acid,
  #4'-ethenyl-2'-hydroxy-1,4',4a-trimethyl-5-oxospiro[2,3,4,7,8,8a-hexahydronaphthalene-6,1'-cyclopentane]-1-carboxylic acid

#####After manual edits, this file was prepared and uploaded to MetaboAnalyst
#####A heatmap was generated in MetaboAnalyst and the image was exported.
#####The image file was used to generate Supplemental Figure 4

##########################End of identification section##########################

##########################End of analysis##########################
