######This analysis generates Figure 3#####

#Summarizing ClassyFire Results

###Read in files after running ClassyFire

#############MSL1
ontMSL1_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesBt_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL1_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesCv_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL1_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesPs_PolarPos_Class.txt",header=FALSE,sep="\t")

ontMSL1_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesBt_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL1_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesCv_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL1_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesPs_PolarNeg_Class.txt",header=FALSE,sep="\t")

#############MSL2
ontMSL2_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL2_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL2_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_PolarPos_Class.txt",header=FALSE,sep="\t")

ontMSL2_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL2_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL2_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_PolarNeg_Class.txt",header=FALSE,sep="\t")

ontMSL2_Bt_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_NonPolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL2_Cv_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_NonPolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL2_Ps_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_NonPolarPos_Class.txt",header=FALSE,sep="\t")

ontMSL2_Bt_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_NonPolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL2_Cv_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_NonPolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL2_Ps_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_NonPolarNeg_Class.txt",header=FALSE,sep="\t")

#############MSL3
ontMSL3_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL3_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_PolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL3_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_PolarPos_Class.txt",header=FALSE,sep="\t")

ontMSL3_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL3_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL3_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg_Class.txt",header=FALSE,sep="\t")

ontMSL3_Bt_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL3_Cv_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos_Class.txt",header=FALSE,sep="\t")
ontMSL3_Ps_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos_Class.txt",header=FALSE,sep="\t")

ontMSL3_Bt_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_NonPolarNeg_Class.txt",header=FALSE,sep="\t")
ontMSL3_Cv_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_NonPolarNeg_Class.txt",header=FALSE,sep="\t")
#ontMSL3_Ps_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_NonPolarNeg_Class.txt",header=FALSE,sep="\t")

BtMSL1 <- rbind(ontMSL1_Bt_PP,ontMSL1_Bt_PN)
BtMSL2 <- rbind(ontMSL2_Bt_PP,ontMSL2_Bt_PN,ontMSL2_Bt_NPP,ontMSL2_Bt_NPN)
BtMSL3 <- rbind(ontMSL3_Bt_PP,ontMSL3_Bt_PN,ontMSL3_Bt_NPP,ontMSL3_Bt_NPN)

Bt1 <- (as.data.frame(table(BtMSL1$V1)))
Bt1$level <- "L1"
Bt2 <- (as.data.frame(table(BtMSL2$V1)))
Bt2$level <- "L2"
Bt3 <- (as.data.frame(table(BtMSL3$V1)))
Bt3$level <- "L3"
BtLevels <- rbind(Bt1,Bt2,Bt3)
BtLevels$org <- "Bt"


CvMSL1 <- rbind(ontMSL1_Cv_PP,ontMSL1_Cv_PN)
CvMSL2 <- rbind(ontMSL2_Cv_PP,ontMSL2_Cv_PN,ontMSL2_Cv_NPP,ontMSL2_Cv_NPN)
CvMSL3 <- rbind(ontMSL3_Cv_PP,ontMSL3_Cv_PN,ontMSL3_Cv_NPP,ontMSL3_Cv_NPN)

Cv1 <- (as.data.frame(table(CvMSL1$V1)))
Cv1$level <- "L1"
Cv2 <- (as.data.frame(table(CvMSL2$V1)))
Cv2$level <- "L2"
Cv3 <- (as.data.frame(table(CvMSL3$V1)))
Cv3$level <- "L3"
CvLevels <- rbind(Cv1,Cv2,Cv3)
CvLevels$org <- "Cv"


PsMSL1 <- rbind(ontMSL1_Ps_PP,ontMSL1_Ps_PN)
PsMSL2 <- rbind(ontMSL2_Ps_PP,ontMSL2_Ps_PN,ontMSL2_Ps_NPP,ontMSL2_Ps_NPN)
PsMSL3 <- rbind(ontMSL3_Ps_PP,ontMSL3_Ps_PN,ontMSL3_Ps_NPP)

Ps1 <- (as.data.frame(table(PsMSL1$V1)))
Ps1$level <- "L1"
Ps2 <- (as.data.frame(table(PsMSL2$V1)))
Ps2$level <- "L2"
Ps3 <- (as.data.frame(table(PsMSL3$V1)))
Ps3$level <- "L3"
PsLevels <- rbind(Ps1,Ps2,Ps3)
PsLevels$org <- "Ps"

AllLevels <- rbind(BtLevels,CvLevels,PsLevels)
#Obtain top 10 categories
library(plyr)
allLevelsSum <- ddply(AllLevels, .(Var1), summarise, Freq=sum(Freq))
library(dplyr)
allLevelsCat <- allLevelsSum %>% top_n(10)

#Obtain top 10 from AllLevels
AllLevelsT10 <- AllLevels[which(AllLevels$Var1 %in% allLevelsCat$Var1),]

#Plot
library(ggplot2)

#Prepare a colorblind palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73")

plot <- ggplot(AllLevelsT10, aes(x=Var1, y=as.numeric(Freq), fill=level)) +
       geom_bar(position=position_dodge2(preserve = 'single'),stat="identity") +
       scale_fill_manual(values=cbPalette)

plots <- plot + facet_grid(rows = vars(org)) + labs(y = "Frequency",fill="Identification Confidence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=16),axis.title.x = element_blank(),axis.text.y = element_text(size=16),axis.title.y = element_text(size=16),legend.title=element_text(size=12))


ggsave("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CLASSYFIRE_ClassLevel_allPolaritiesAndIonsCombined.eps",plot=plots,device="eps",width=30, units="cm",dpi=600)

###Direct Parent level

#############MSL1
ontMSL1_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesBt_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL1_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesCv_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL1_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesPs_PolarPos_DirectParent.txt",header=FALSE,sep="\t")

ontMSL1_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesBt_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL1_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesCv_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL1_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_JGIrefMatchesPs_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")

#############MSL2
ontMSL2_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_PolarPos_DirectParent.txt",header=FALSE,sep="\t")

ontMSL2_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")

ontMSL2_Bt_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Cv_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Ps_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")

ontMSL2_Bt_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesBt_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Cv_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesCv_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL2_Ps_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_MonaNistMatchesPs_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")

#############MSL3
ontMSL3_Bt_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Cv_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_PolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Ps_PP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_PolarPos_DirectParent.txt",header=FALSE,sep="\t")

ontMSL3_Bt_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Cv_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Ps_PN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_PolarNeg_DirectParent.txt",header=FALSE,sep="\t")

ontMSL3_Bt_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Cv_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Ps_NPP <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_NonPolarPos_DirectParent.txt",header=FALSE,sep="\t")

ontMSL3_Bt_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesBt_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")
ontMSL3_Cv_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesCv_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")
#ontMSL3_Ps_NPN <- read.csv("MassSpec/ClassyFire/initialFiles/ClassyFireOutput/ClassyFire_SIRIUSFinIDMatchesPs_NonPolarNeg_DirectParent.txt",header=FALSE,sep="\t")

BtMSL1 <- rbind(ontMSL1_Bt_PP,ontMSL1_Bt_PN)
BtMSL2 <- rbind(ontMSL2_Bt_PP,ontMSL2_Bt_PN,ontMSL2_Bt_NPP,ontMSL2_Bt_NPN)
BtMSL3 <- rbind(ontMSL3_Bt_PP,ontMSL3_Bt_PN,ontMSL3_Bt_NPP,ontMSL3_Bt_NPN)

Bt1 <- (as.data.frame(table(BtMSL1$V1)))
Bt1$level <- "L1"
Bt2 <- (as.data.frame(table(BtMSL2$V1)))
Bt2$level <- "L2"
Bt3 <- (as.data.frame(table(BtMSL3$V1)))
Bt3$level <- "L3"
BtLevels <- rbind(Bt1,Bt2,Bt3)
BtLevels$org <- "Bt"


CvMSL1 <- rbind(ontMSL1_Cv_PP,ontMSL1_Cv_PN)
CvMSL2 <- rbind(ontMSL2_Cv_PP,ontMSL2_Cv_PN,ontMSL2_Cv_NPP,ontMSL2_Cv_NPN)
CvMSL3 <- rbind(ontMSL3_Cv_PP,ontMSL3_Cv_PN,ontMSL3_Cv_NPP,ontMSL3_Cv_NPN)

Cv1 <- (as.data.frame(table(CvMSL1$V1)))
Cv1$level <- "L1"
Cv2 <- (as.data.frame(table(CvMSL2$V1)))
Cv2$level <- "L2"
Cv3 <- (as.data.frame(table(CvMSL3$V1)))
Cv3$level <- "L3"
CvLevels <- rbind(Cv1,Cv2,Cv3)
CvLevels$org <- "Cv"


PsMSL1 <- rbind(ontMSL1_Ps_PP,ontMSL1_Ps_PN)
PsMSL2 <- rbind(ontMSL2_Ps_PP,ontMSL2_Ps_PN,ontMSL2_Ps_NPP,ontMSL2_Ps_NPN)
PsMSL3 <- rbind(ontMSL3_Ps_PP,ontMSL3_Ps_PN,ontMSL3_Ps_NPP)

Ps1 <- (as.data.frame(table(PsMSL1$V1)))
Ps1$level <- "L1"
Ps2 <- (as.data.frame(table(PsMSL2$V1)))
Ps2$level <- "L2"
Ps3 <- (as.data.frame(table(PsMSL3$V1)))
Ps3$level <- "L3"
PsLevels <- rbind(Ps1,Ps2,Ps3)
PsLevels$org <- "Ps"

AllLevels <- rbind(BtLevels,CvLevels,PsLevels)
#Obtain top 10 categories
library(plyr)
allLevelsSum <- ddply(AllLevels, .(Var1), summarise, Freq=sum(Freq))
library(dplyr)
allLevelsCat <- allLevelsSum %>% top_n(10)

#Obtain top 10 from AllLevels
AllLevelsT10 <- AllLevels[which(AllLevels$Var1 %in% allLevelsCat$Var1),]

#Plot
library(ggplot2)
plot <- ggplot(AllLevelsT10, aes(x=Var1, y=as.numeric(Freq), fill=level)) +
       geom_bar(position=position_dodge2(preserve = 'single'),stat="identity")
plots <- plot + facet_grid(rows = vars(org)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/MonoCulturePaper/MassSpec/ClassyFire/CLASSYFIRE_DirectParentLevel_allPolaritiesAndIonsCombined.tif",plot=plots,device="tiff",width=30, units="cm",dpi=600)
