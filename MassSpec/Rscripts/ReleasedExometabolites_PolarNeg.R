#####This is analysis is determine the number of released exometabolites in each isolate#####
#####and, determine the identity of a subset of these exometabolites#####
#####This is for Polar analysis in negative ionization mode#####

library(dplyr)

#Read in MZmine final feature table for MS analysis
#####################MAKE SURE YOU CHANGE THIS ON GITHUB FOR FINAL CODE SUBMISSION#####################
filt <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegwExConForQC.csv",check.names=F)

#Remove outlier-1mem_Bt_45hr-4
filt_f <- filt[,-c(84)]

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
write.csv(exoMetabMonoCulture,"MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",row.names=FALSE)
exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)
#Order columns
exoMetabolomeMonoCultureOrdered <- exoMetabolomeMonoCulture[ , order(names(exoMetabolomeMonoCulture))]
write.csv(exoMetabolomeMonoCultureOrdered,"MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonocultureOrdered.csv",row.names=FALSE)

#####Start of released exometabolite analysis#####
exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)

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
#18611 features remain

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
groupedMonoSum$conditions <- rep(c("Bt","Cv","Ps"),each=6,18611)
#Add time
groupedMonoSum$time <- rep(c(12.5,25,30,35,40,45),55833) #This number is the number above * 3 (# of conditions)

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

BtMinMaxIDs <- unique(BtmonoFC$ID[BtmonoFC$FC>3]) #2136
CvMinMaxIDs <- unique(CvmonoFC$ID[CvmonoFC$FC>3]) #1147
PsMinMaxIDs <- unique(PsmonoFC$ID[PsmonoFC$FC>3]) #1802

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

#Remove features whose lowest RSD value is still >0.2. High variable features
BtFinalIDs <- groupRSDSumBtMAX$ID[groupRSDSumBtMAX$RSD<0.2] #1327
CvFinalIDs <- groupRSDSumCvMAX$ID[groupRSDSumCvMAX$RSD<0.2] #638
PsFinalIDs <- groupRSDSumPsMAX$ID[groupRSDSumPsMAX$RSD<0.2] #975

#Confirm filtered IDs that also had a RSD value of 0 was the RSD value at the initial timepoint
groupRSDSumBT0s[which(groupRSDSumBT0s$ID %in% BtFinalIDs),] #Remove ID 50612
groupRSDSumCv0s[which(groupRSDSumCv0s$ID %in% CvFinalIDs),] #Looks good
groupRSDSumPs0s[which(groupRSDSumPs0s$ID %in% PsFinalIDs),] #One feature in 25hr but also at initial TP, too. This is okay

BtFinalIDs <- BtFinalIDs[!BtFinalIDs %in% 50612] #1326

#Save IDs for each isolate that have been determined to be released
saveRDS(BtFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
saveRDS(CvFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
saveRDS(PsFinalIDs, "MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

##########End of filtering steps to determine released exometabolites##########





###################Start here to continue with identification###################
BtFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
CvFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
PsFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

#Create a venn diagram
finalFeats <- list(Bt=BtFinalIDsPolarNeg,Cv=CvFinalIDsPolarNeg,Ps=PsFinalIDsPolarNeg)
finalFeatsV <- lapply(finalFeats, function(x) x[!is.na(x)])
library(VennDiagram)

#Supplementary Figure 1
venn.diagram(finalFeatsV,filename="Figures/SFig1/MonoOverlappingFeatures_PolarNeg.tiff")

#Venn output shows
#1050 unique features for Bt
#372 unique features for Cv
#662 unique features for Ps

#####Use venn diagrams from Polar/NonPolar and Positive/Negative modes (Supplemental Figure 1) to generate Table 1#####

#Let's export a CSV file for all released exometabolites for Metaboanalyst

#Read in file
exoMetabolome <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsPolarNeg,CvFinalIDsPolarNeg,PsFinalIDsPolarNeg))
#Remove NA
finalFeatsJGI <- finalFeatsJGI[!is.na(finalFeatsJGI)]

#Create CSV for all secreted exometabolites
exoMetabolomeReleased <- exoMetabolome[c(which(exoMetabolome$ID %in% unique(finalFeatsJGI))),]
#Add internal standard- 13C-15N-Proline, ID 363
exoMetabolomeReleased <- rbind(exoMetabolomeReleased,exoMetabolome[exoMetabolome$ID=="363",])

#Re-order columns
exoMetabolomeReleasedO <- exoMetabolomeReleased[ , order(names(exoMetabolomeReleased))]
#Write csv
write.csv(exoMetabolomeReleasedO,"MassSpec/releaseAnalysis/MS/outputFiles/MZminePolarNegallReleasedMetabolitesIndividualSamples.csv")

#Note: MZminePolarNegallReleasedMetabolitesIndividualSamples.csv file was then manually edited to prepare for Metaboanalyst.
#This was done is two ways:
#1) MZminePolarNegallReleasedMetabolitesIndividualSamples_readyForMetaboAnalyst_wTime.csv for HeatMap analysis.
    #Perform statistical analysis on MetaboAnalyst and generate heatmap averagign bioreps: https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml
    #Note, pseudosamples were added to the data where only two bioreps were present. This was performed because MetaboAnalyst requries at least 3 samples per group.
    #The pseudosample was simply an average of the two bioreps, so it has no effect on the overall results after averaging bioreps
#2) MZminePolarNegallReleasedMetabolitesIndividualSamples_readyForMetaboAnalyst.csv for PCA analysis
    #Data was normalized in metaboanalyst and normalized data was exported to PolarNeg_allReleased_IndSamples_MetaboanalystNormalization.csv
      #This was then uploaded in R to perform PCA in vegan, generating Figure 2.

#MZmineNonPolarNegallReleasedMetabolitesIndividualSamples.csv manually edited file was uploade to MetaboAnalyst
#This file was used to generate Figure 1 and Supplemental Figure 2
#Follow the Rscript, "PCAs.R", to finish generating Figure 1.
#Supp. Fig. 2 was generated directly in MetaboAnalyst





#####################Identification of released metabolites#####################
BtFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
CvFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
PsFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)

finalFeatsJGI <- unique(c(BtFinalIDsPolarNeg,CvFinalIDsPolarNeg,PsFinalIDsPolarNeg))

filtJGI <- exoMetabolomeMonoCulture[c(which(exoMetabolomeMonoCulture$ID %in% unique(finalFeatsJGI))),]

###Steps for exometabolite identification
# 1) JGI Ref library
# 2) Nist17
# 3) MoNa
#For any overlaps, #1 overrides. For overlaps between 2 & 3, highest cosine similarity overrides

###########Identification step 1) JGI Ref library###########

#Upload JGI library
library(dplyr)
JGIrefLib <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/JGI_refLibrary_PolarNegUpdated.csv",check.names=F)

#Extract released metabolites with JGI MSL1 identification

allJGI=NULL
for(i in 1:nrow(JGIrefLib)){
  mzLow= JGIrefLib$mzMin[i]
  mzHigh= JGIrefLib$mzMax[i]
  rtLow= JGIrefLib$rtMin[i]
  rtHigh= JGIrefLib$rtMax[i]

  new_frame<- exoMetabolomeMonoCulture %>% filter(between(MZ,mzLow,mzHigh)) %>% filter(between(RT,rtLow,rtHigh))
  print(nrow(new_frame))
  if(nrow(new_frame) > 0){
  if(nrow(new_frame) > 1){
  RTdev <- abs(new_frame$RT-JGIrefLib$avgRt[i])
  print(RTdev)
  ID <- new_frame$ID[which(RTdev == min(RTdev))]
  Label <- JGIrefLib$Label[i]
  df <- data.frame(ID,Label)
  allJGI=rbind(allJGI,df)
  } else{
  ID <- new_frame$ID[1]
  Label <- JGIrefLib$Label[i]
  df <- data.frame(ID,Label)
  allJGI=rbind(allJGI,df)
  }
  } else {
  }
}

#Obtain JGI IDed compounds that pass criteria for released
allJGIreleased <- allJGI[which(allJGI$ID %in% filtJGI$ID),]

#Check which ID had isomers
table(table(allJGI$Label)) #32 0's, 141 1's- 32 exometabolites not detected
table(table(allJGI$ID))  #Of the 141 identified, 103 unique, 19 isomers (103 + 2*19 =141)
#Extract all JGI MSL1 exometabolites
allMSL1 <- exoMetabolomeMonoCulture[which(exoMetabolomeMonoCulture$ID %in% allJGI$ID),]
#Order columns
allMSL1.O <- allMSL1[ , order(names(allMSL1))]
write.csv(allMSL1.O,"MassSpec/releaseAnalysis/MSMS/outputFiles/PolarNeg_JGI_MSL1_IDs.csv",row.names=FALSE)
#Export label names
write.csv(allJGI,"MassSpec/releaseAnalysis/MSMS/outputFiles/PolarNeg_JGI_MSL1_Labels.csv",row.names=FALSE)

#####Use PolarPos_JGI_MSL1_IDs.csv to find the top most accumulated exometabolites#####

#Extract released Bt exometabolites from JGI
Btmetabs <- allJGIreleased[which(allJGIreleased$ID %in% BtFinalIDsPolarNeg),]

#Extract released Cv exometabolites from JGI
Cvmetabs <- allJGIreleased[which(allJGIreleased$ID %in% CvFinalIDsPolarNeg),]

#Extract released Ps exometabolites from JGI
Psmetabs <- allJGIreleased[which(allJGIreleased$ID %in% PsFinalIDsPolarNeg),]

###########Identification step 2) Nist17 & MoNa###########

#Now, let's match these IDs with those that have MS/MS data or matches to JGI standard library. Having an MS/MS does not mean it was an identified JGI MSL1 match.

#Load MSMS analysis from MZmine
MSMSdat <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZmine_PolarNegMSMS_peakArea.csv",check.names=F)

#They match! Filtered features with MSMS data
BtMSMS <- BtFinalIDsPolarNeg[which(BtFinalIDsPolarNeg %in% MSMSdat$"row ID")] #235/1326 filtered features
CvMSMS <- CvFinalIDsPolarNeg[which(CvFinalIDsPolarNeg %in% MSMSdat$"row ID")] #193/638 filtered features
PsMSMS <- PsFinalIDsPolarNeg[which(PsFinalIDsPolarNeg %in% MSMSdat$"row ID")] #205/975 filtered features

#See which IDs match with MSMS data
MSMSmatchesNist17 <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZminePolarNegMSMSnist17Matches.csv",check.names=F)
MSMSmatchesMONA <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/PolarNegMSMS_MONAspectraMatches.csv",check.names=F)

#MSMS matches from nist17
BtMSSpecsNist17 <- BtMSMS[which(BtMSMS %in% MSMSmatchesNist17$"cluster index")]
CvMSSpecsNist17 <- CvMSMS[which(CvMSMS %in% MSMSmatchesNist17$"cluster index")]
PsMSSpecsNist17 <- PsMSMS[which(PsMSMS %in% MSMSmatchesNist17$"cluster index")]

MSSpecsNist17 <- unique(c(BtMSSpecsNist17, CvMSSpecsNist17,PsMSSpecsNist17))
#Remove those found from step 1
n17remove <- MSSpecsNist17[which(MSSpecsNist17 %in% allJGIreleased$ID)]
MSSpecsNist17F <- MSSpecsNist17[!MSSpecsNist17 %in% n17remove]
which(MSSpecsNist17F %in% allJGIreleased$ID) #nothing. Good.

#MSMS matches from MONA
BtMSSpecsMONA <- BtMSMS[which(BtMSMS %in% MSMSmatchesMONA$rowID)]
CvMSSpecsMONA <- CvMSMS[which(CvMSMS %in% MSMSmatchesMONA$rowID)]
PsMSSpecsMONA <- PsMSMS[which(PsMSMS %in% MSMSmatchesMONA$rowID)]

MSSpecsMONA <- unique(c(BtMSSpecsMONA, CvMSSpecsMONA,PsMSSpecsMONA))
#No matches at all!

#Remove those found from step 1
#monAremove <- MSSpecsMONA[which(MSSpecsMONA %in% allJGIreleased$ID)]
#MSSpecsMONAF <- MSSpecsMONA[!MSSpecsMONA %in% monAremove]
#which(MSSpecsMONAF %in% allJGIreleased$ID) #nothing. Good.

#Of those identified by nist17 or MONA, which ones overlap?
#MvN <- MSSpecsNist17F[which(MSSpecsNist17F %in% MSSpecsMONAF)]

#N17 <- MSMSmatchesNist17[which(MSMSmatchesNist17$"cluster index" %in% MvN),]
#MONA <- MSMSmatchesMONA[which(MSMSmatchesMONA$rowID %in% MvN),]

#Order rows
#N17o <- N17[order(N17[,11]),]
#MONAo <- MONA[with(MONA,order(rowID)),]

#Which identifications to use from Nist17?
#monaRemove <- N17o[N17o$MQScore>MONAo$cos,11]
#Which identifications to use from MONA?
#nist17Remove <- MONAo[MONAo$cos>N17o$MQScore,1]

#Remove NIST17 and MONA overlaps
#MSSpecsMONAFf <- MSSpecsMONAF[!MSSpecsMONAF %in% monaRemove]
#MSSpecsNist17Ff <- MSSpecsNist17F[!MSSpecsNist17F %in% nist17Remove]
#which(MSSpecsNist17Ff %in% MSSpecsMONAFf) #none. Good.

#No overlaps found

#MONAhits <- MSMSmatchesMONA[which(MSMSmatchesMONA$rowID %in% MSSpecsMONAF),]
#library(stringr)
#MONAhits$LibraryID <- str_split_fixed(MONAhits$molecule, "as", 2)[,1]
nist17hits <- MSMSmatchesNist17[which(MSMSmatchesNist17$"cluster index" %in% MSSpecsNist17F),]

#Now we need to rid of redundancies within a method
#For some reason, the mzMine procedure picks redundant peaks. We want to keep the ID with the higher cosine score

#Check for redundant compounds

#JGI redundants
length(allJGIreleased$ID) #69
length(unique(allJGIreleased$ID)) #58
length(allJGIreleased$Label) #69
length(unique(allJGIreleased$Label)) #69
#Note, redundancies are typically isomers.
allJGIreleased$ID[which(duplicated(allJGIreleased$ID))]

#Nist17 redundants
length(nist17hits$"cluster index")   #13
length(unique(nist17hits$"cluster index"))  #13
length(nist17hits$LibraryID)   #13
length(unique(nist17hits$LibraryID))  #12
#obtain redundant features
dupeNist17 <- nist17hits[,"LibraryID"]
nist17hits[duplicated(dupeNist17) | duplicated(dupeNist17, fromLast=TRUE),]

#MONA redundants
#length(MONAhits$rowID) #0
#length(unique(MONAhits$rowID)) #0
#length(MONAhits$LibraryID) #0
#length(unique(MONAhits$LibraryID)) #0

#Subset columns
#MONAhitsSub <- MONAhits[,c(1,6)]
nist17hitsSub <- nist17hits[,c(10,2)]

#Add additional column for source of identification
#MONAhitsSub$source <- rep("mona",nrow(MONAhitsSub))
nist17hitsSub$source <- rep("nist17",nrow(nist17hitsSub))
allJGIreleased$source <- rep("JGI",nrow(allJGIreleased))

#Rename columns
#names(MONAhitsSub) <- c("ID","Label","source")
names(nist17hitsSub) <- c("ID","Label","source")

#Bring it all together
identifiedFinal <- rbind(allJGIreleased,nist17hitsSub)

#Obtain IDS from all methods
filtFinalMSMS <- exoMetabolomeMonoCulture[c(which(exoMetabolomeMonoCulture$ID %in% unique(identifiedFinal$ID))),]
#Add internal standard- 13C-15N-Proline, ID 363
filtFinalMSMS <- rbind(filtFinalMSMS,exoMetabolomeMonoCulture[exoMetabolomeMonoCulture$ID=="363",])
#Re-order columns
filtFinalMSMSO <- filtFinalMSMS[ , order(names(filtFinalMSMS))]
#Write csv for individual samples
write.csv(filtFinalMSMSO,"MassSpec/releaseAnalysis/MSMS/outputFiles/MZminePolarNegReleasedMetaboliteswIDS_IndividualSamples.csv")

#Write compounds to CSV
write.csv(identifiedFinal,"MassSpec/releaseAnalysis/MSMS/outputFiles/PolarNegMSMScompounds.csv")

######################################################Manually edit file######################################################

#The file: MZminePolarNegReleasedMetaboliteswIDS_IndividualSamples.csv was manually edited

#First, align and copy over identified metabolites and source from PolarNegMSMScompounds.csv onto the MZminePolarNegReleasedMetaboliteswIDS_IndividualSamples CSV file
#Remember, JGI IDs are isomers
allJGIreleased$ID[which(duplicated(allJGIreleased$ID))]

#remove similar annotated metabolites from Mona and Nist17 against the JGI reference
toMatch <- as.character(JGIrefLib$Label)
#To get an idea of the matches
#matchesMONA <- unique(grep(paste(toMatch,collapse="|"), MONAhitsSub$Label ,ignore.case = T, value=TRUE))
matchesNIST17 <- unique(grep(paste(toMatch,collapse="|"), nist17hitsSub$Label ,ignore.case = T, value=TRUE))

###NOTE, this is a rough match. Do not get rid of all these "matches". For example, 2'-deoxyadenosine shows for matchesMONA
###because JGI has 2'-DEOXYGUANOSINE. 2'-deoxy is the match. Don't get rid of 2'-deoxyadenosine

#####And, manually remove redunancies within the same identification #Keep the exometabolite with the higher cosine score
#MONAhits[duplicated(dupeMONA) | duplicated(dupeMONA, fromLast=TRUE),]
nist17hits[duplicated(dupeNist17) | duplicated(dupeNist17, fromLast=TRUE),]

#Lastly, look through molecules. Is it logical that it should be there? Specifically for MONA and NIST17 hits.
#For example, remove 3,7,8,2'-Tetrahydroxyflavone,Testosterone, Threonic acid,  ###FINISH THIS INSPECTION

#####After manual edits, this file was prepared and uploaded to MetaboAnalyst
#####A heatmap was generated in MetaboAnalyst and the image was exported.
#####The image file was used to generate Supplemental Figure 4

##########################End of identification section##########################





##################Finding most accumulated MSL1 exometabolites##################

JGIPN_IDs <- read.csv("MassSpec/releaseAnalysis/MSMS/outputFiles/PolarNeg_JGI_MSL1_IDs.csv",sep=",",header=TRUE,check.names=FALSE)
exoMetabolomeMonoCulture <- read.csv("MassSpec/releaseAnalysis/MS/initialFiles/MZminePolarNegMonoculture.csv",check.names=F)

#Avg BioReps
library(reshape2)
JGIPN_IDsMelt <- melt(JGIPN_IDs,id.vars=c("ID","RT","MZ"))

#Change NA values and 0s to low value
JGIPN_IDsMelt[,5][JGIPN_IDsMelt[,5]<2] <- 2

#We need to split names so we can average conditions within a time point

JGIsplit <- data.frame(do.call('rbind', strsplit(as.character(JGIPN_IDsMelt$variable),'-',fixed=TRUE)))
#Place this into melted object
JGIPN_IDsMelt$groupTime <- JGIsplit$X1

library(dplyr)
detach("package:plyr", unload=TRUE)
JGIPN_IDsMelt_grouped <- JGIPN_IDsMelt %>% group_by(ID,RT,MZ,groupTime)
#Get mean
JGIPN_IDsMelt_groupedAvg <- summarise(JGIPN_IDsMelt_grouped, mean=mean(value))

#Convert back to DF
library(reshape2)
JGIPN_Avg.df <- dcast(JGIPN_IDsMelt_groupedAvg, ID + MZ + RT ~ groupTime)

#Get top 5 most abundant exometabolites for B.thailandensis at 45 hr
Bt_exoMetab <- JGIPN_Avg.df[,c(1:3,9)]
Bt_mostAbunO <- Bt_exoMetab[order(-(Bt_exoMetab$"1mem_Bt_45hr")),]
#Obtain IDs of top five
Bt_mostAbunTF <-Bt_mostAbunO$ID[1:5]
#Bt contains a redundant exometabolite (4-pyridoxate) in both pos and neg mode.
#The next most abundant metabolite in pos or neg mode happens to be in negative mode. We'll swap this ID with the third ID (the redundant 4-pyridoxate)
Bt_mostAbunTF <-Bt_mostAbunO$ID[c(1:2,4:5,6)]

#Get top 5 most abundant exometabolites for C.violaceum at 45 hr
Cv_exoMetab <- JGIPN_Avg.df[,c(1:3,15)]
Cv_mostAbunO <- Cv_exoMetab[order(-(Cv_exoMetab$"1mem_Cv_45hr")),]
#Obtain IDs of top five
Cv_mostAbunTF <-Cv_mostAbunO$ID[1:5]

#Get top 5 most abundant exometabolites for P.syringae at 45 hr
Ps_exoMetab <- JGIPN_Avg.df[,c(1:3,21)]
Ps_mostAbunO <- Ps_exoMetab[order(-(Ps_exoMetab$"1mem_Ps_45hr")),]
#Obtain IDs of top five
Ps_mostAbunTF <-Ps_mostAbunO$ID[1:5]

#Match IDs to JGI compounds

#Upload JGI library
library(dplyr)
JGIrefLib <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/JGI_refLibrary_PolarNegUpdated.csv",check.names=F)

#Extract all MS1 metabolites from JGI

allJGI=NULL
for(i in 1:nrow(JGIrefLib)){
  mzLow= JGIrefLib$mzMin[i]
  mzHigh= JGIrefLib$mzMax[i]
  rtLow= JGIrefLib$rtMin[i]
  rtHigh= JGIrefLib$rtMax[i]

  new_frame<- exoMetabolomeMonoCulture %>% filter(between(MZ,mzLow,mzHigh)) %>% filter(between(RT,rtLow,rtHigh))
  print(nrow(new_frame))
  if(nrow(new_frame) > 0){
  if(nrow(new_frame) > 1){
  RTdev <- abs(new_frame$RT-JGIrefLib$avgRt[i])
  print(RTdev)
  ID <- new_frame$ID[which(RTdev == min(RTdev))]
  Label <- JGIrefLib$Label[i]
  df <- data.frame(ID,Label)
  allJGI=rbind(allJGI,df)
  } else{
  ID <- new_frame$ID[1]
  Label <- JGIrefLib$Label[i]
  df <- data.frame(ID,Label)
  allJGI=rbind(allJGI,df)
  }
  } else {
  }
}

Bt_PN_topFiveMSL1 <- allJGI$Label[which(allJGI$ID %in% Bt_mostAbunTF)]
Cv_PN_topFiveMSL1 <- allJGI$Label[which(allJGI$ID %in% Cv_mostAbunTF)]
Ps_PN_topFiveMSL1 <- allJGI$Label[which(allJGI$ID %in% Ps_mostAbunTF)]

#Of these top accumulated exometabolites, which fit the criteria of released?
BtFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
CvFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
PsFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

BtmostAbunAndReleased <- allJGI$Label[which(allJGI$ID %in% BtFinalIDsPolarNeg[which(BtFinalIDsPolarNeg %in% Bt_mostAbunTF)])]
CvmostAbunAndReleased <- allJGI$Label[which(allJGI$ID %in% CvFinalIDsPolarNeg[which(CvFinalIDsPolarNeg %in% Cv_mostAbunTF)])]
PsmostAbunAndReleased <- allJGI$Label[which(allJGI$ID %in% PsFinalIDsPolarNeg[which(PsFinalIDsPolarNeg %in% Ps_mostAbunTF)])]

#This was used to generate Figure 4

##########################End of analysis##########################
