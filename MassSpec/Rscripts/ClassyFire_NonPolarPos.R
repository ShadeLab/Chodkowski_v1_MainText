#####This analysis is to determine the chemical ontology of released exometabolites in each isolate#####
#####This is for Nonpolar analysis in positive ionization mode#####

#####CSIFingerID Prep#####
BtFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsNonPolarPos.rds")
CvFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsNonPolarPos.rds")
PsFinalIDsNonPolarPos <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsNonPolarPos.rds")

MSMSdat <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZmine_NonPolarPosMSMS_peakArea.csv",check.names=F)

#They match! Filtered features with MSMS data
BtMSMS <- BtFinalIDsNonPolarPos[which(BtFinalIDsNonPolarPos %in% MSMSdat$"row ID")] #153/383 filtered features
CvMSMS <- CvFinalIDsNonPolarPos[which(CvFinalIDsNonPolarPos %in% MSMSdat$"row ID")] #147/551 filtered features
PsMSMS <- PsFinalIDsNonPolarPos[which(PsFinalIDsNonPolarPos %in% MSMSdat$"row ID")] #92/623 filtered features

#Read in the final table of IDed metabolites. This has been manually manipulated to remove any additional redundancies between MSL2 data
finalMetaboliteNPP <- read.csv("MassSpec/releaseAnalysis/MSMS/outputFiles/manualEdits/MZmineNonPolarPosReleasedMetaboliteswIDS_IndividualSamples_manualEdits.csv", header=TRUE,sep=",",check.names = FALSE)

#One last Venn for identified metabolites
BtMSSpecsIDed <- BtMSMS[which(BtMSMS %in% finalMetaboliteNPP$ID)]
CvMSSpecsIDed <- CvMSMS[which(CvMSMS %in% finalMetaboliteNPP$ID)]
PsMSSpecsIDed <- PsMSMS[which(PsMSMS %in% finalMetaboliteNPP$ID)]

#After manually curating the final set of secreted metabolites, let's go back and find out which came from which organism
#Extract Bt secreted metabolites
Btmetabs <- finalMetaboliteNPP[which(finalMetaboliteNPP$ID %in% BtFinalIDsNonPolarPos),]

#Extract Cv secreted metabolites
Cvmetabs <- finalMetaboliteNPP[which(finalMetaboliteNPP$ID %in% CvFinalIDsNonPolarPos),]

#Extract Ps secreted metabolites
Psmetabs <- finalMetaboliteNPP[which(finalMetaboliteNPP$ID %in% PsFinalIDsNonPolarPos),]

###Extract MSL2 data###
monaNist <- read.csv("MassSpec/ClassyFire/initialFiles/MonaNist_MSL2_forCSIID_NonPolarPos.csv",header=TRUE,sep=",")

#Filter MONA/Nist only
MSMSmatchesMONABt <- Btmetabs[Btmetabs$source=="mona",]
MSMSmatchesMONACv <- Cvmetabs[Cvmetabs$source=="mona",]
MSMSmatchesMONAPs <- Psmetabs[Psmetabs$source=="mona",]

MSMSmatchesNISTBt <- Btmetabs[Btmetabs$source=="nist17",]
MSMSmatchesNISTCv <- Cvmetabs[Cvmetabs$source=="nist17",]
MSMSmatchesNISTPs <- Psmetabs[Psmetabs$source=="nist17",]

#Make sure all labels match across
match(MSMSmatchesMONABt$ID,monaNist$ID) #good
match(MSMSmatchesMONACv$ID,monaNist$ID) #good
match(MSMSmatchesMONAPs$ID,monaNist$ID) #good

match(MSMSmatchesNISTBt$ID,monaNist$ID) #good
match(MSMSmatchesNISTCv$ID,monaNist$ID) #good
match(MSMSmatchesNISTPs$ID,monaNist$ID) #good

#Obtain matches
MONAMatchesBt <- monaNist[which(monaNist$ID %in% MSMSmatchesMONABt$ID),]
MONAMatchesCv <- monaNist[which(monaNist$ID %in% MSMSmatchesMONACv$ID),]
MONAMatchesPs <- monaNist[which(monaNist$ID %in% MSMSmatchesMONAPs$ID),]

NISTMatchesBt <- monaNist[which(monaNist$ID %in% MSMSmatchesNISTBt$ID),]
NISTMatchesCv <- monaNist[which(monaNist$ID %in% MSMSmatchesNISTCv$ID),]
NISTMatchesPs <- monaNist[which(monaNist$ID %in% MSMSmatchesNISTPs$ID),]

#Combine
msL2Bt <- rbind(MONAMatchesBt,NISTMatchesBt)
msL2Cv <- rbind(MONAMatchesCv,NISTMatchesCv)
msL2Ps <- rbind(MONAMatchesPs,NISTMatchesPs)

#Remove dups- ready to export.
msL2BtF <- msL2Bt[!duplicated(msL2Bt$ID),]
msL2CvF <- msL2Cv[!duplicated(msL2Cv$ID),]
msL2PsF <- msL2Ps[!duplicated(msL2Ps$ID),]





###########Skip this section. The CSIFingerIDResults results are too large to upload to GitHub#########

#####################Data is available upon request##################

###Extract MSL3###
#Obtain secreted metabolites that were not identified at MSL1/MSL2
BtmetabsNOID <- BtFinalIDsNonPolarPos[!BtFinalIDsNonPolarPos %in% Btmetabs$ID]
CvmetabsNOID <- CvFinalIDsNonPolarPos[!CvFinalIDsNonPolarPos %in% Cvmetabs$ID]
PsmetabsNOID <- PsFinalIDsNonPolarPos[!PsFinalIDsNonPolarPos %in% Psmetabs$ID]

#Obtain these from the SIRIUS:CSIFingerID
for(i in 1:length(BtmetabsNOID)){
fileC <- BtmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/BtNonPolarPosCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/BtNonPolarPosSIRIUSFINID_MSL3")
}
else{
}
}

for(i in 1:length(CvmetabsNOID)){
fileC <- CvmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/CvNonPolarPosCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/CvNonPolarPosSIRIUSFINID_MSL3")
}
else{
}
}

for(i in 1:length(PsmetabsNOID)){
fileC <- PsmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/PsNonPolarPosCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/PsNonPolarPosSIRIUSFINID_MSL3")
}
else{
}
}

#Export for MSL2
write.csv(msL2BtF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/NonPolarPos_ForClassyFire_msL2Bt.csv")
write.csv(msL2CvF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/NonPolarPos_ForClassyFire_msL2Cv.csv")
write.csv(msL2PsF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/NonPolarPos_ForClassyFire_msL2Ps.csv")
#Export for MSL3
  #In shell, Extract the second line (the top CSIFingerID hit) from each file.

sed -s -n 2p MassSpec/ClassyFire/outputFiles/BtNonPolarPosSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/BtNonPolarPosSIRIUSFINID_MSL3.txt
sed -s -n 2p MassSpec/ClassyFire/outputFiles/CvNonPolarPosSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/CvNonPolarPosSIRIUSFINID_MSL3.txt
sed -s -n 2p MassSpec/ClassyFire/outputFiles/PsNonPolarPosSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PsNonPolarPosSIRIUSFINID_MSL3.txt

##############################MANUAL EDITS##############################################################

#Need to run classyfire: http://classyfire.wishartlab.com/
#For each MS classification level for each organism.
#After running classyfire, download the SDF file.

#Re-name files to correspond to sample of interest. (E.g. 27367.sdf > ClassyFire_JGIrefMatchesBt_NonPolarPos.sdf)
#Extract class and direct parent information

awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesBt_NonPolarPos.sdf > ClassyFire_MonaNistMatchesBt_NonPolarPos_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesCv_NonPolarPos.sdf > ClassyFire_MonaNistMatchesCv_NonPolarPos_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesPs_NonPolarPos.sdf > ClassyFire_MonaNistMatchesPs_NonPolarPos_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos_Class.txt

awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesBt_NonPolarPos.sdf > ClassyFire_MonaNistMatchesBt_NonPolarPos_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesCv_NonPolarPos.sdf > ClassyFire_MonaNistMatchesCv_NonPolarPos_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesPs_NonPolarPos.sdf > ClassyFire_MonaNistMatchesPs_NonPolarPos_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos.sdf > ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos_DirectParent.txt

#Place txt files in directory of interest

#################################End of manual edits#################################################################

###########Done with skipped section. Analysis complete#########

###########Continue with "PlottingClassyFireResults" script to create Figure 3#########
