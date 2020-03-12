#####This analysis is for differential gene expression######
#####And for creating transcript cumulative sum curves######
#####Finally, transporters for each isolate are analyzed for differential expression and expression above a low threshold#####

#############################Burkholderia##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Bt_diffExp_all.csv",header=TRUE,sep=",")
mappingMono <- mapping[mapping$class=="con",]

#remove non-numerical columns from count matrix
genecounts <- counts[5:ncol(counts)]
row.names(genecounts) <- counts[,1]

#Obtain only mono samples
libMono <- lib[lib$Condition=="Mono",]

genecountsMono <- genecounts[,which(colnames(genecounts) %in% libMono$libraryName)]

#Replace NA values with 0s
genecountsMono[is.na(genecountsMono)] <- 0

#Remove 0 count genes
genesZeroC <- rownames(genecountsMono)[which(rowSums(genecountsMono) == 0)]
genesRemove <- c(genesZeroC)

#Extract unique identifiers from each condition
libMono$libraryName = as.character(libMono$libraryName)
genecounts.sort=genecountsMono[,match(libMono$libraryName, names(genecountsMono))]
library(stringr)
Extract <- c("mono")
keywords <- str_extract(libMono$sampleName, paste(Extract,collapse="|"))

all <- which(keywords %in% c("mono"))
samps <- as.vector(libMono$sampleName[all])
samples <- match(samps,libMono$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mappingMono[,1])

#Align sample names between objects
diffGenes =geneMat[,match(align, colnames(geneMat))]
colnames(diffGenes) <- mappingMono$Code
rownames(diffGenes) <- rownames(genecounts.sort)

#Remove genes of interest
diffGenes = diffGenes[-c(which(rownames(diffGenes) %in% genesRemove)),]
#5641 > 5633 genes

#Now, let's remove genes with less than 10 counts in 90% of the samples- 94 samples so 85 is ~90% of samples

idx <- rowSums(diffGenes >= 10 ) >= 20


diffGenesF <- diffGenes[idx,]
#5634 > 5175 genes

#Prepare files for DESEq analysis
mappingMono$time <- as.factor(as.numeric(as.character(mappingMono$time)))
ddsMat <- DESeqDataSetFromMatrix(countData = diffGenesF, colData=mappingMono, design= ~time)
ddsMat$time <- relevel(ddsMat$time, ref = "750")
#ddsMat$time <- factor(ddsMat$time, levels=c("750","1500","1800","2100","2400","2700"))
ddsMat <- DESeq(ddsMat,betaPrior=FALSE)
saveRDS(ddsMat, file = "RNAseq/diffGeneExp/outputFiles/BtmonoDeseqObject.rds")
BtmonoDeseq <- results(ddsMat)
saveRDS(BtmonoDeseq, file = "RNAseq/diffGeneExp/outputFiles/BtmonoDeseq.rds")

############End OF differential gene expression analysis for B.thailandensis############





############Create transcript cumulative sum curve for B.thailandensis to determine low expression threshold############

selectNorm <- counts(ddsMat,normalized=TRUE)

#Extract 45 TP
mapping45 <- mappingMono[mappingMono$time=="2700",]
selectNorm45 <- selectNorm[,which(colnames(selectNorm) %in% mapping45$Code)]

selectNorm45DF <- as.data.frame(rowMeans(selectNorm45))
colnames(selectNorm45DF) <- c("mean")
selectNorm45DF$geneID <- row.names(selectNorm45DF)

#orderedHtoL <- selectNorm45DF[order(-selectNorm45DF$mean)]
orderedHtoL <- selectNorm45DF[order(-selectNorm45DF$mean),]
orderedHtoL$transcript <- rep(1,nrow(orderedHtoL))

plot(cumsum(orderedHtoL$transcript)/sum(orderedHtoL$transcript),cumsum(orderedHtoL$mean)/sum(orderedHtoL$mean))

#orderedHtoLow <- orderedHtoL[3106:nrow(orderedHtoL),]
#orderedHightoL <- orderedHtoL[1:157,]
#orderedHtoLowPercent <- sum(orderedHtoLow$mean)/sum(orderedHtoL$mean)
#quantile(orderedHtoL$mean[order(orderedHtoL$mean)], c(.03,.60))

quantile(orderedHtoL$mean, c(.25,.97))

###Quantile 25%- 158.9032

#############################Chromobacterium##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cv_diffExp_all.csv",header=TRUE,sep=",")
mappingMono <- mapping[mapping$class=="con",]

#remove non-numerical columns from count matrix
genecounts <- counts[5:ncol(counts)]
row.names(genecounts) <- counts[,1]

#Obtain only mono samples
libMono <- lib[lib$Condition=="Mono",]

genecountsMono <- genecounts[,which(colnames(genecounts) %in% libMono$libraryName)]

#Replace NA values with 0s
genecountsMono[is.na(genecountsMono)] <- 0

#Remove 0 count genes
genesZeroC <- rownames(genecountsMono)[which(rowSums(genecountsMono) == 0)]
genesRemove <- c(genesZeroC)

#Extract unique identifiers from each condition
libMono$libraryName = as.character(libMono$libraryName)
genecounts.sort=genecountsMono[,match(libMono$libraryName, names(genecountsMono))]
library(stringr)
Extract <- c("mono")
keywords <- str_extract(libMono$sampleName, paste(Extract,collapse="|"))

all <- which(keywords %in% c("mono"))
samps <- as.vector(libMono$sampleName[all])
samples <- match(samps,libMono$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mappingMono[,1])

#Align sample names between objects
diffGenes =geneMat[,match(align, colnames(geneMat))]
colnames(diffGenes) <- mappingMono$Code
rownames(diffGenes) <- rownames(genecounts.sort)

#Remove genes of interest
diffGenes = diffGenes[-c(which(rownames(diffGenes) %in% genesRemove)),]
#4371 > 4369 genes

#Now, let's remove genes with less than 10 counts in 90% of the samples- 94 samples so 85 is ~90% of samples

idx <- rowSums(diffGenes >= 10 ) >= 21


diffGenesF <- diffGenes[idx,]
#4369> 3805 genes

#Prepare files for DESEq analysis
mappingMono$time <- as.factor(as.numeric(as.character(mappingMono$time)))
#mappingMono$rep <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,6,1,2,3,4,5,6)
#mappingMono$rep <- as.factor(as.numeric(as.character(mappingMono$rep)))

ddsMat <- DESeqDataSetFromMatrix(countData = diffGenesF, colData=mappingMono, design= ~time)
ddsMat$time <- relevel(ddsMat$time, ref = "750")
#ddsMat$time <- factor(ddsMat$time, levels=c("750","1500","1800","2100","2400","2700"))
ddsMat <- DESeq(ddsMat, betaPrior=FALSE)
saveRDS(ddsMat, file = "RNAseq/diffGeneExp/outputFiles/CvmonoDeseqObject.rds")
CvmonoDeseq <- results(ddsMat)
saveRDS(CvmonoDeseq, file = "RNAseq/diffGeneExp/outputFiles/CvmonoDeseq.rds")

############End OF differential gene expression analysis for C.violaceum############





############Create transcript cumulative sum curve for C.violaceum to determine low expression threshold############

selectNorm <- counts(ddsMat,normalized=TRUE)

#Extract 45 TP
mapping45 <- mappingMono[mappingMono$time=="2700",]
selectNorm45 <- selectNorm[,which(colnames(selectNorm) %in% mapping45$Code)]

selectNorm45DF <- as.data.frame(rowMeans(selectNorm45))
colnames(selectNorm45DF) <- c("mean")
selectNorm45DF$geneID <- row.names(selectNorm45DF)

#orderedHtoL <- selectNorm45DF[order(-selectNorm45DF$mean)]
orderedHtoL <- selectNorm45DF[order(-selectNorm45DF$mean),]
orderedHtoL$transcript <- rep(1,nrow(orderedHtoL))

plot(cumsum(orderedHtoL$transcript)/sum(orderedHtoL$transcript),cumsum(orderedHtoL$mean)/sum(orderedHtoL$mean))

orderedHtoLow <- orderedHtoL[3106:nrow(orderedHtoL),]

#orderedHightoL <- orderedHtoL[1:157,]
#orderedHtoLowPercent <- sum(orderedHtoLow$mean)/sum(orderedHtoL$mean)
#orderedHightoLPercent <-
#quantile(orderedHtoL$mean[order(orderedHtoL$mean)], c(.03,.60))

quantile(orderedHtoL$mean, c(.25,.97))

#25% quantile- 142.768

#############################Pseudomonas##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Psraw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Ps_diffExp_all.csv",header=TRUE,sep=",")
mappingMono <- mapping[mapping$class=="con",]

#remove non-numerical columns from count matrix
genecounts <- counts[5:ncol(counts)]
row.names(genecounts) <- counts[,1]

#Obtain only mono samples
libMono <- lib[lib$Condition=="Mono",]

genecountsMono <- genecounts[,which(colnames(genecounts) %in% libMono$libraryName)]

#Replace NA values with 0s
genecountsMono[is.na(genecountsMono)] <- 0

#Remove 0 count genes
genesZeroC <- rownames(genecountsMono)[which(rowSums(genecountsMono) == 0)]
genesRemove <- c(genesZeroC)

#Extract unique identifiers from each condition
libMono$libraryName = as.character(libMono$libraryName)
genecounts.sort=genecountsMono[,match(libMono$libraryName, names(genecountsMono))]
library(stringr)
Extract <- c("mono")
keywords <- str_extract(libMono$sampleName, paste(Extract,collapse="|"))

all <- which(keywords %in% c("mono"))
samps <- as.vector(libMono$sampleName[all])
samples <- match(samps,libMono$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mappingMono[,1])

#Align sample names between objects
diffGenes =geneMat[,match(align, colnames(geneMat))]
colnames(diffGenes) <- mappingMono$Code
rownames(diffGenes) <- rownames(genecounts.sort)

#Remove genes of interest
diffGenes = diffGenes[-c(which(rownames(diffGenes) %in% genesRemove)),]
#5853 > 5835 genes

#Now, let's remove genes with less than 10 counts in 90% of the samples- 94 samples so 85 is ~90% of samples

idx <- rowSums(diffGenes >= 10 ) >= 22


diffGenesF <- diffGenes[idx,]
#5835> 5536 genes

#Prepare files for DESEq analysis
mappingMono$time <- as.factor(as.numeric(as.character(mappingMono$time)))
ddsMat <- DESeqDataSetFromMatrix(countData = diffGenesF, colData=mappingMono, design= ~time)
ddsMat$time <- relevel(ddsMat$time, ref = "750")
#ddsMat$time <- factor(ddsMat$time, levels=c("750","1500","1800","2100","2400","2700"))
ddsMat <- DESeq(ddsMat,betaPrior=FALSE)
saveRDS(ddsMat, file = "RNAseq/diffGeneExp/outputFiles/PsmonoDeseqObject.rds")
PsmonoDeseq <- results(ddsMat)
saveRDS(PsmonoDeseq, file = "RNAseq/diffGeneExp/outputFiles/PsmonoDeseq.rds")

############End OF differential gene expression analysis for P.syringae############





############Create transcript cumulative sum curve for P.syringae to determine low expression threshold############

selectNorm <- counts(ddsMat,normalized=TRUE)

#Extract 45 TP
mapping45 <- mappingMono[mappingMono$time=="2700",]
selectNorm45 <- selectNorm[,which(colnames(selectNorm) %in% mapping45$Code)]

selectNorm45DF <- as.data.frame(rowMeans(selectNorm45))
colnames(selectNorm45DF) <- c("mean")
selectNorm45DF$geneID <- row.names(selectNorm45DF)

orderedHtoL <- selectNorm45DF[order(-selectNorm45DF$mean),]
orderedHtoL$transcript <- rep(1,nrow(orderedHtoL))

plot(cumsum(orderedHtoL$transcript)/sum(orderedHtoL$transcript),cumsum(orderedHtoL$mean)/sum(orderedHtoL$mean))

#orderedHtoLow <- orderedHtoL[3106:nrow(orderedHtoL),]
#orderedHightoL <- orderedHtoL[1:157,]
#orderedHtoLowPercent <- sum(orderedHtoLow$mean)/sum(orderedHtoL$mean)
#orderedHightoLPercent <-
#quantile(orderedHtoL$mean[order(orderedHtoL$mean)], c(.03,.60))

quantile(orderedHtoL$mean, c(.25,.97))

#25% quantile- 174.2475


################Transport analysis- Is expression above the 25% quantile threshold at TP 45?#########################################

#############################Burkholderia##############################
library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Bt_diffExp_all.csv",header=TRUE,sep=",")

transPorters <- read.csv("RNAseq/diffGeneExp/initialFiles/BthailandensisE264TransAAPTransporters.csv",header=TRUE,sep=",")

#load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses/lfc/output/BthailandensisLFC.RData")
countsTransP <- counts[which(counts$geneID %in% transPorters$GeneID),]

countsTransP45h <- countsTransP[,which(colnames(countsTransP) %in% c("BXNYU","BXOAA","BXOBN"))]
countsTransP45h$mean <- rowMeans(countsTransP45h)
countsTransP45h$GeneID <- countsTransP$geneID

#From above, the 25% quantile for the cumulative abundance curve at TP45 for Burkholderia:
#25% quantile- 158.9032

#How many genes are above this threshold?
table(countsTransP45h$mean>158.9032)
countsTransP45hFiltered <- countsTransP45h[countsTransP45h$mean>158.9032,]

#Extract these genes from TransportFile
BttransportLowEThresh <- transPorters[which(transPorters$GeneID %in% countsTransP45hFiltered$GeneID),]
write.csv(BttransportLowEThresh,"RNAseq/diffGeneExp/outputFiles/BttransportLowEThresh.csv")

BtmonoDeseq <- readRDS("RNAseq/diffGeneExp/outputFiles/BtmonoDeseq.rds")
#obtain fdr values less than 0.01
table(BtmonoDeseq$padj<0.01)
BtmonoFDR <- BtmonoDeseq$padj<0.01

#Convert NAs to FALSE
BtmonoFDR[is.na(BtmonoFDR)]<- FALSE

#Obtain rows with pval less than 0.01
BtmonoFDRcutOff <- BtmonoDeseq[BtmonoFDR,]

#Which transport genes are significantly diff regulated?
BtFinalTransportGenes <- countsTransP45hFiltered[which(countsTransP45hFiltered$GeneID %in% rownames(BtmonoFDRcutOff)),]
BtFinalTransportGenesID <- transPorters[which(transPorters$GeneID %in% BtFinalTransportGenes$GeneID),]
write.csv(BtFinalTransportGenesID,"RNAseq/diffGeneExp/outputFiles/BttransportLowEThreshAndpADJ.csv")


#############################Chromobacterium##############################
library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cv_diffExp_all.csv",header=TRUE,sep=",")

transPorters <- read.csv("RNAseq/diffGeneExp/initialFiles/Cviolaceum31532TransAAPTransporters.csv",header=TRUE,sep=",")

#load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses/lfc/output/CviolaceumLFC.RData")
countsTransP <- counts[which(counts$geneID %in% transPorters$GeneID),]

countsTransP45h <- countsTransP[,which(colnames(countsTransP) %in% c("BZTUU","BZTXA","BZTYN","BZTZU"))]
countsTransP45h$mean <- rowMeans(countsTransP45h)
countsTransP45h$GeneID <- countsTransP$geneID

#From above, the 25% quantile for the cumulative abundance curve at TP45 for Chromobacterium:
#25% quantile- 142.768

#How many genes are above this threshold?
table(countsTransP45h$mean>142.768)
countsTransP45hFiltered <- countsTransP45h[countsTransP45h$mean>142.768,]

#Extract these genes from TransportFile
CvtransportLowEThresh <- transPorters[which(transPorters$GeneID %in% countsTransP45hFiltered$GeneID),]
write.csv(CvtransportLowEThresh,"RNAseq/diffGeneExp/outputFiles/CvtransportLowEThresh.csv")

CvmonoDeseq <- readRDS("RNAseq/diffGeneExp/outputFiles/CvmonoDeseq.rds")
#obtain fdr values less than 0.01
table(CvmonoDeseq$padj<0.01)
CvmonoFDR <- CvmonoDeseq$padj<0.01

#Convert NAs to FALSE
CvmonoFDR[is.na(CvmonoFDR)]<- FALSE

#Obtain rows with pval less than 0.01
CvmonoFDRcutOff <- CvmonoDeseq[CvmonoFDR,]

#Which transport genes are significantly diff regulated?
CvFinalTransportGenes <- countsTransP45hFiltered[which(countsTransP45hFiltered$GeneID %in% rownames(CvmonoFDRcutOff)),]
CvFinalTransportGenesID <- transPorters[which(transPorters$GeneID %in% CvFinalTransportGenes$GeneID),]
write.csv(CvFinalTransportGenesID,"RNAseq/diffGeneExp/outputFiles/CvtransportLowEThreshAndpADJ.csv")

#############################Pseudomonas##############################
library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Psraw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Ps_diffExp_all.csv",header=TRUE,sep=",")

transPorters <- read.csv("RNAseq/diffGeneExp/initialFiles/PsyringaeDC3000TransAAPTransporters.csv",header=TRUE,sep=",")

#load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses/lfc/output/PsyringaeLFC.RData")
countsTransP <- counts[which(counts$geneID %in% transPorters$GeneID),]

countsTransP45h <- countsTransP[,which(colnames(countsTransP) %in% c("BYSXY","BYSZG","BYTAS","BYTBY"))]
countsTransP45h$mean <- rowMeans(countsTransP45h)
countsTransP45h$GeneID <- countsTransP$geneID

#From above, the 25% quantile for the cumulative abundance curve at TP45 for Chromobacterium:
#25% quantile- 174.2475

#How many genes are above this threshold?
table(countsTransP45h$mean>174.2475)
countsTransP45hFiltered <- countsTransP45h[countsTransP45h$mean>174.2475,]

#Extract these genes from TransportFile
PstransportLowEThresh <- transPorters[which(transPorters$GeneID %in% countsTransP45hFiltered$GeneID),]
write.csv(PstransportLowEThresh,"RNAseq/diffGeneExp/outputFiles/PstransportLowEThresh.csv")

PsmonoDeseq <- readRDS("RNAseq/diffGeneExp/outputFiles/PsmonoDeseq.rds")
#obtain fdr values less than 0.01
table(PsmonoDeseq$padj<0.01)
PsmonoFDR <- PsmonoDeseq$padj<0.01

#Convert NAs to FALSE
PsmonoFDR[is.na(PsmonoFDR)]<- FALSE

#Obtain rows with pval less than 0.01
PsmonoFDRcutOff <- PsmonoDeseq[PsmonoFDR,]

#Which transport genes are significantly diff regulated?
PsFinalTransportGenes <- countsTransP45hFiltered[which(countsTransP45hFiltered$GeneID %in% rownames(PsmonoFDRcutOff)),]
PsFinalTransportGenesID <- transPorters[which(transPorters$GeneID %in% PsFinalTransportGenes$GeneID),]
write.csv(PsFinalTransportGenesID,"RNAseq/diffGeneExp/outputFiles/PstransportLowEThreshAndpADJ.csv")

######From the annotated transporters for each isolate, an excel sheet containing loci in involved in transport were color-coded as follows:

#Isoltate_transportLowEThreshAndpADJ.csv used to assess criteria for being expression above low threshold and differentially expressed- If passed, colored blue
#Isolate_transportLowEThresh.csv used to assess criteria for being expression above low threshold and differentially expressed- If passed, colored yellow
#Loci not passing either criteria were colored red

#Final excel sheets were placed in RNAseq/diffGeneExp/outputFiles/manualEdits

###############################Supplemental File X ###############################
