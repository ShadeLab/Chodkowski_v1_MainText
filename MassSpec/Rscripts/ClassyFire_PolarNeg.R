#####This analysis is to determine the chemical ontology of released exometabolites in each isolate#####
#####This is for Polar analysis in negative ionization mode#####

#####CSIFingerID Prep#####
BtFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/BtFinalIDsPolarNeg.rds")
CvFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/CvFinalIDsPolarNeg.rds")
PsFinalIDsPolarNeg <- readRDS("MassSpec/releaseAnalysis/MS/outputFiles/PsFinalIDsPolarNeg.rds")

MSMSdat <- read.csv("MassSpec/releaseAnalysis/MSMS/initialFiles/MZmine_PolarNegMSMS_peakArea.csv",check.names=F)

#They match! Filtered features with MSMS data
BtMSMS <- BtFinalIDsPolarNeg[which(BtFinalIDsPolarNeg %in% MSMSdat$"row ID")] #235/1327 filtered features
CvMSMS <- CvFinalIDsPolarNeg[which(CvFinalIDsPolarNeg %in% MSMSdat$"row ID")] #193/638 filtered features
PsMSMS <- PsFinalIDsPolarNeg[which(PsFinalIDsPolarNeg %in% MSMSdat$"row ID")] #205/975 filtered features

#Read in the final table of IDed metabolites. This has been manually manipulated to remove any additional redundancies between MSL1 data MSL2
finalMetabolitePN <- read.csv("MassSpec/releaseAnalysis/MSMS/outputFiles/manualEdits/MZminePolarNegReleasedMetaboliteswIDS_IndividualSamples_manualEdits.csv", header=TRUE,sep=",",check.names = FALSE)

#One last Venn for only identified metabolites
BtMSSpecsIDed <- BtMSMS[which(BtMSMS %in% finalMetabolitePN$ID)]
CvMSSpecsIDed <- CvMSMS[which(CvMSMS %in% finalMetabolitePN$ID)]
PsMSSpecsIDed <- PsMSMS[which(PsMSMS %in% finalMetabolitePN$ID)]

#After manually curating the final set of secreted metabolites, let's go back and find out which came from which organism
#Extract Bt secreted metabolites
Btmetabs <- finalMetabolitePN[which(finalMetabolitePN$ID %in% BtFinalIDsPolarNeg),]

#Extract Cv secreted metabolites
Cvmetabs <- finalMetabolitePN[which(finalMetabolitePN$ID %in% CvFinalIDsPolarNeg),]

#Extract Ps secreted metabolites
Psmetabs <- finalMetabolitePN[which(finalMetabolitePN$ID %in% PsFinalIDsPolarNeg),]

###Extract MSL1 data###
JGIref <- read.csv("MassSpec/ClassyFire/initialFiles/JGIRefLib_MSL1_forCSIID_PolarNeg.csv", header=TRUE,sep=",")

#Filter JGI only
MSMSmatchesJGIBt <- Btmetabs[Btmetabs$source=="JGI",]
MSMSmatchesJGICv <- Cvmetabs[Cvmetabs$source=="JGI",]
MSMSmatchesJGIPs <- Psmetabs[Psmetabs$source=="JGI",]

#Make sure all labels match across
match(MSMSmatchesJGIBt$Label,JGIref$Label) #good
match(MSMSmatchesJGICv$Label,JGIref$Label) #good
match(MSMSmatchesJGIPs$Label,JGIref$Label) #good

#Obtain matches- ready to export
JGIrefMatchesBt <- JGIref[which(JGIref$Label %in% MSMSmatchesJGIBt$Label),]
JGIrefMatchesCv <- JGIref[which(JGIref$Label %in% MSMSmatchesJGICv$Label),]
JGIrefMatchesPs <- JGIref[which(JGIref$Label %in% MSMSmatchesJGIPs$Label),]

###Extract MSL2 data###
monaNist <- read.csv("MassSpec/ClassyFire/initialFiles/MonaNist_MSL2_forCSIID_PolarNeg.csv",header=TRUE,sep=",")

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
BtmetabsNOID <- BtFinalIDsPolarNeg[!BtFinalIDsPolarNeg %in% Btmetabs$ID]
CvmetabsNOID <- CvFinalIDsPolarNeg[!CvFinalIDsPolarNeg %in% Cvmetabs$ID]
PsmetabsNOID <- PsFinalIDsPolarNeg[!PsFinalIDsPolarNeg %in% Psmetabs$ID]

#Obtain these from the SIRIUS:CSIFingerID
for(i in 1:length(BtmetabsNOID)){
fileC <- BtmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/BtPolarNegCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/BtPolarNegSIRIUSFINID_MSL3")
}
else{
}
}

for(i in 1:length(CvmetabsNOID)){
fileC <- CvmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/CvPolarNegCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/CvPolarNegSIRIUSFINID_MSL3")
}
else{
}
}

for(i in 1:length(PsmetabsNOID)){
fileC <- PsmetabsNOID[i]
checkFile <- paste0("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CSIFingerIDResults/PsPolarNegCSIFinID/",fileC,".csv")
if(file.exists(checkFile) == TRUE){
file.copy(checkFile,"MassSpec/ClassyFire/outputFiles/PsPolarNegSIRIUSFINID_MSL3")
}
else{
}
}

#Export for MSL1
write.csv(JGIrefMatchesBt,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_JGIrefMatchesBt.csv")
write.csv(JGIrefMatchesCv,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_JGIrefMatchesCv.csv")
write.csv(JGIrefMatchesPs,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_JGIrefMatchesPs.csv")
#Export for MSL2
write.csv(msL2BtF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_msL2Bt.csv")
write.csv(msL2CvF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_msL2Cv.csv")
write.csv(msL2PsF,"MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PolarNeg_ForClassyFire_msL2Ps.csv")
#Export for MSL3
  #In shell, Extract the second line (the top CSIFingerID hit) from each file.

sed -s -n 2p MassSpec/ClassyFire/outputFiles/BtPolarNegSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/BtPolarNegSIRIUSFINID_MSL3.txt
sed -s -n 2p MassSpec/ClassyFire/outputFiles/CvPolarNegSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/CvPolarNegSIRIUSFINID_MSL3.txt
sed -s -n 2p MassSpec/ClassyFire/outputFiles/PsPolarNegSIRIUSFINID_MSL3/*.csv > MassSpec/ClassyFire/outputFiles/forClassyFireAnalysis/PsPolarNegSIRIUSFINID_MSL3.txt

##############################MANUAL EDITS##############################################################

#Need to run classyfire: http://classyfire.wishartlab.com/
#For each MS classification level for each organism.
#After running classyfire, download the SDF file.

#Re-name files to correspond to sample of interest. (E.g. 27367.sdf > ClassyFire_JGIrefMatchesBt_PolarNeg.sdf)
#Extract class and direct parent information

awk 'f{print;f=0} /Class/{f=1}' ClassyFire_JGIrefMatchesBt_PolarNeg.sdf > ClassyFire_JGIrefMatchesBt_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_JGIrefMatchesCv_PolarNeg.sdf > ClassyFire_JGIrefMatchesCv_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_JGIrefMatchesPs_PolarNeg.sdf > ClassyFire_JGIrefMatchesPs_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesBt_PolarNeg.sdf > ClassyFire_MonaNistMatchesBt_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesCv_PolarNeg.sdf > ClassyFire_MonaNistMatchesCv_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_MonaNistMatchesPs_PolarNeg.sdf > ClassyFire_MonaNistMatchesPs_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg_Class.txt
awk 'f{print;f=0} /Class/{f=1}' ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg_Class.txt

awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_JGIrefMatchesBt_PolarNeg.sdf > ClassyFire_JGIrefMatchesBt_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_JGIrefMatchesCv_PolarNeg.sdf > ClassyFire_JGIrefMatchesCv_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_JGIrefMatchesPs_PolarNeg.sdf > ClassyFire_JGIrefMatchesPs_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesBt_PolarNeg.sdf > ClassyFire_MonaNistMatchesBt_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesCv_PolarNeg.sdf > ClassyFire_MonaNistMatchesCv_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_MonaNistMatchesPs_PolarNeg.sdf > ClassyFire_MonaNistMatchesPs_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg_DirectParent.txt
awk 'f{print;f=0} /Direct Parent/{f=1}' ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg.sdf > ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg_DirectParent.txt

#Place txt files in directory of interest

#################################End of manual edits#################################################################

###########Done with skipped section. Analysis complete#########

###########Continue with "PlottingClassyFireResults" script to create Figure 3#########
