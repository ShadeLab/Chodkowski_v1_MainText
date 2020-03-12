library(DESeq2)
library(pathview)
library(gage)
library(KEGGREST)


##################################Burkholderia##################################

lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Psraw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Ps_diffExp_all.csv",header=TRUE,sep=",")

#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/PsmonoDeseqObject.rds")

#LFC data
load("RNAseq/LFC/outputFiles/BthailandensisLFC.RData")

norm <- counts(ddsMat, normalized=TRUE)
geneID <- counts[which(counts$geneID %in% row.names(norm)),]

paths <- unique(keggLink("pathway","bte"))
pathsv <- gsub("\\D+","", paths)
#Remove global and overview maps and, pathways with less than 3 nodes
global <- c("01100","01110","01120","01200","01210","01212","01230","01220","00511","01054")
keggPaths <- pathsv[-c(match(global,pathsv))]

###Monoculture

lfc_mono_25.fc <- lfc_mono_25withinComparisionToTime0.fc
lfc_mono_30.fc <- lfc_mono_30withinComparisionToTime0.fc
lfc_mono_35.fc <- lfc_mono_35withinComparisionToTime0.fc
lfc_mono_40.fc <- lfc_mono_40withinComparisionToTime0.fc
lfc_mono_45.fc <- lfc_mono_45withinComparisionToTime0.fc

lfc_mono_df <- data.frame(cbind(lfc_mono_25.fc,lfc_mono_30.fc,lfc_mono_35.fc,lfc_mono_40.fc,lfc_mono_45.fc))
row.names(lfc_mono_df) <- geneID$Locus.Tag
lfc_mono.df <- as.matrix(lfc_mono_df)

#obtain lfc on kegg pathway maps for monoculture

setwd("RNAseq/KeggPathway/Bthailandensis")
pv.out <- pathview(gene.data = lfc_mono.df[, 1:5], pathway.id = keggPaths,species = "bte", out.suffix = keggPaths,multi.state=T,gene.idtype="kegg")



##################################Chromobacterium##################################

library(DESeq2)
library(pathview)
library(gage)
library(KEGGREST)

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cv_diffExp_all.csv",header=TRUE,sep=",")

KI <- read.csv(file="RNAseq/KeggPathway/Cviolaceum_blastKoala_TheFinal.csv",sep=",",header=TRUE)

#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/CvmonoDeseqObject.rds")

#LFC data
load("RNAseq/LFC/outputFiles/CviolaceumLFC.RData")

norm <- counts(ddsMat, normalized=TRUE)
geneID <- counts[which(counts$geneID %in% row.names(norm)),]

paths <- unique(keggLink("pathway","cvi"))
pathsv <- gsub("\\D+","", paths)
#Remove global and overview maps and, pathways with less than 3 nodes orelse the program stalls
#Pathways with 3 nodes or less- 00511,
global <- c("01100","01110","01120","01200","01210","01212","01230","01220","00511","01054")
keggPaths <- pathsv[-c(match(global,pathsv))]


###Monoculture

lfc_mono_25.fc <- lfc_mono25_withinComparision_ToTime0$log2FoldChange
lfc_mono_30.fc <- lfc_mono30_withinComparision_ToTime0$log2FoldChange
lfc_mono_35.fc <- lfc_mono35_withinComparision_ToTime0$log2FoldChange
lfc_mono_40.fc <- lfc_mono40_withinComparision_ToTime0$log2FoldChange
lfc_mono_45.fc <- lfc_mono45_withinComparision_ToTime0$log2FoldChange

lfc_mono_df <- data.frame(cbind(lfc_mono_25.fc,lfc_mono_30.fc,lfc_mono_35.fc,lfc_mono_40.fc,lfc_mono_45.fc))
row.names(lfc_mono_df) <- geneID$geneID

#Now we need to convert geneID to KEGG IDS
lfc_mono_df <- lfc_mono_df[which(row.names(lfc_mono_df) %in% KI$Query),]
row.names(lfc_mono_df) <- KI$KO[match(row.names(lfc_mono_df),KI$Query)]

lfc_mono.df <- as.matrix(lfc_mono_df)

#obtain lfc on kegg pathway maps for monoculture

setwd("RNAseq/KeggPathway/Cviolaceum")
pv.out <- pathview(gene.data = lfc_mono.df[, 1:5], pathway.id = keggPaths,species = "ko", out.suffix = keggPaths,multi.state=T,gene.idtype="kegg",kegg.native=T)



##################################Pseudomonas##################################

library(DESeq2)
library(pathview)
library(gage)
library(KEGGREST)

lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Psraw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Ps_diffExp_all.csv",header=TRUE,sep=",")

KI <- read.csv(file="RNAseq/KeggPathway/Ps_blastKoala_TheFinal.csv",sep=",",header=TRUE)

#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/PsmonoDeseqObject.rds")

#LFC data
load("RNAseq/LFC/outputFiles/PsyringaeLFC.RData")

norm <- counts(ddsMat, normalized=TRUE)
geneID <- counts[which(counts$geneID %in% row.names(norm)),]

paths <- unique(keggLink("pathway","pst"))
pathsv <- gsub("\\D+","", paths)
#Remove global and overview maps and, pathways with less than 3 nodes
global <- c("01100","01110","01120","01200","01210","01212","01230","01220")
keggPaths <- pathsv[-c(match(global,pathsv))]


###Monoculture

lfc_mono_25.fc <- lfc_mono25_withinComparision_ToTime0$log2FoldChange
lfc_mono_30.fc <- lfc_mono30_withinComparision_ToTime0$log2FoldChange
lfc_mono_35.fc <- lfc_mono35_withinComparision_ToTime0$log2FoldChange
lfc_mono_40.fc <- lfc_mono40_withinComparision_ToTime0$log2FoldChange
lfc_mono_45.fc <- lfc_mono45_withinComparision_ToTime0$log2FoldChange

lfc_mono_df <- data.frame(cbind(lfc_mono_25.fc,lfc_mono_30.fc,lfc_mono_35.fc,lfc_mono_40.fc,lfc_mono_45.fc))
row.names(lfc_mono_df) <- geneID$geneID

#Now we need to convert geneID to KEGG IDS
lfc_mono_df <- lfc_mono_df[which(row.names(lfc_mono_df) %in% KI$Query),]
row.names(lfc_mono_df) <- KI$KO[match(row.names(lfc_mono_df),KI$Query)]

lfc_mono.df <- as.matrix(lfc_mono_df)

#obtain lfc on kegg pathway maps for monoculture

setwd("RNAseq/KeggPathway/Psyringae")
pv.out <- pathview(gene.data = lfc_mono.df[, 1:5], pathway.id = keggPaths,species = "ko", out.suffix = keggPaths,multi.state=T,gene.idtype="kegg",kegg.native=T)

################End of KEGG Pathway analysis##########################

#####These data were used to create pathways involved in succinate production#####
