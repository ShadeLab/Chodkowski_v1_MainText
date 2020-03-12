#####
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS: /opt/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/lib/libR.so
LAPACK: /opt/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/modules/lapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] DESeq2_1.22.2               SummarizedExperiment_1.12.0
 [3] DelayedArray_0.8.0          BiocParallel_1.16.1
 [5] matrixStats_0.54.0          Biobase_2.42.0
 [7] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1
 [9] IRanges_2.16.0              S4Vectors_0.20.1
[11] BiocGenerics_0.28.0         ashr_2.2-7

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.5.0          foreach_1.4.4
 [4] Formula_1.2-3          assertthat_0.2.0       latticeExtra_0.6-28
 [7] blob_1.1.1             GenomeInfoDbData_1.2.0 RSQLite_2.1.1
[10] pillar_1.3.1           backports_1.1.2        lattice_0.20-35
[13] glue_1.3.0             digest_0.6.18          RColorBrewer_1.1-2
[16] XVector_0.22.0         checkmate_1.8.5        colorspace_1.3-2
[19] htmltools_0.3.6        Matrix_1.2-14          plyr_1.8.4
[22] XML_3.98-1.16          pkgconfig_2.0.2        genefilter_1.64.0
[25] zlibbioc_1.28.0        purrr_0.2.5            xtable_1.8-3
[28] scales_1.0.0           htmlTable_1.12         tibble_1.4.2
[31] annotate_1.60.0        ggplot2_3.1.0          nnet_7.3-12
[34] lazyeval_0.2.1         survival_2.42-3        magrittr_1.5
[37] crayon_1.3.4           memoise_1.1.0          doParallel_1.0.14
[40] MASS_7.3-50            foreign_0.8-70         truncnorm_1.0-8
[43] tools_3.5.0            data.table_1.12.0      stringr_1.3.1
[46] locfit_1.5-9.1         munsell_0.5.0          cluster_2.0.7-1
[49] AnnotationDbi_1.44.0   bindrcpp_0.2.2         compiler_3.5.0
[52] rlang_0.3.1            grid_3.5.0             RCurl_1.95-4.11
[55] iterators_1.0.10       rstudioapi_0.9.0       htmlwidgets_1.3
[58] bitops_1.0-6           base64enc_0.1-3        gtable_0.2.0
[61] codetools_0.2-15       DBI_1.0.0              R6_2.3.0
[64] gridExtra_2.3          knitr_1.20             dplyr_0.7.8
[67] bit_1.1-14             bindr_0.1.1            Hmisc_4.1-1
[70] stringi_1.2.4          pscl_1.5.2             SQUAREM_2017.10-1
[73] Rcpp_1.0.0             geneplotter_1.60.0     rpart_4.1-13
[76] acepack_1.4.1          tidyselect_0.2.5
#####

#############################Burkholderia##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/B-thailandensis_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Bt_diffExp_all.csv",header=TRUE,sep=",")
#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/BtmonoDeseqObject.rds")

###Obtain log fold changes for a within a condition comparing to the first time point.

#Mono fold changes within
lfc_mono25_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1500","750"))
lfc_mono30_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1800","750"))
lfc_mono35_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2100","750"))
lfc_mono40_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2400","750"))
lfc_mono45_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2700","750"))


#Bt monoculture extraction of lfc at  each time point compared to t initial
lfc_mono_25withinComparisionToTime0.fc <- lfc_mono25_withinComparision_ToTime0$log2FoldChange
lfc_mono_30withinComparisionToTime0.fc <- lfc_mono30_withinComparision_ToTime0$log2FoldChange
lfc_mono_35withinComparisionToTime0.fc <- lfc_mono35_withinComparision_ToTime0$log2FoldChange
lfc_mono_40withinComparisionToTime0.fc <- lfc_mono40_withinComparision_ToTime0$log2FoldChange
lfc_mono_45withinComparisionToTime0.fc <- lfc_mono45_withinComparision_ToTime0$log2FoldChange

#set working directory to place of interest
save.image("RNAseq/LFC/outputFiles/BthailandensisLFC.RData")



#############################Chromobacterium##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cviolaceum_raw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Cv_diffExp_all.csv",header=TRUE,sep=",")
#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/CvmonoDeseqObject.rds")


###Obtain log fold changes for a within a condition comparing to the first time point.

#Mono fold changes within
lfc_mono25_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1500","750"))
lfc_mono30_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1800","750"))
lfc_mono35_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2100","750"))
lfc_mono40_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2400","750"))
lfc_mono45_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2700","750"))


#Cv monoculture extraction of lfc at  each time point compared to t initial
lfc_mono_25withinComparisionToTime0.fc <- lfc_mono25_withinComparision_ToTime0$log2FoldChange
lfc_mono_30withinComparisionToTime0.fc <- lfc_mono30_withinComparision_ToTime0$log2FoldChange
lfc_mono_35withinComparisionToTime0.fc <- lfc_mono35_withinComparision_ToTime0$log2FoldChange
lfc_mono_40withinComparisionToTime0.fc <- lfc_mono40_withinComparision_ToTime0$log2FoldChange
lfc_mono_45withinComparisionToTime0.fc <- lfc_mono45_withinComparision_ToTime0$log2FoldChange

save.image("RNAseq/LFC/outputFiles/CviolaceumLFC.RData")



#############################Pseudomonas##############################

library(DESeq2)
lib <- read.csv(file="RNAseq/diffGeneExp/initialFiles/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Psraw_counts.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="RNAseq/diffGeneExp/initialFiles/Ps_diffExp_all.csv",header=TRUE,sep=",")
#Load ddsMat OBJECT
ddsMat <- readRDS("RNAseq/diffGeneExp/outputFiles/PsmonoDeseqObject.rds")

###Obtain log fold changes for a within a condition comparing to the first time point.

#Mono fold changes within
lfc_mono25_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1500","750"))
lfc_mono30_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "1800","750"))
lfc_mono35_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2100","750"))
lfc_mono40_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2400","750"))
lfc_mono45_withinComparision_ToTime0 <- results(ddsMat,contrast=c("time", "2700","750"))


#Ps monoculture extraction of lfc at  each time point compared to t initial
lfc_mono_25withinComparisionToTime0.fc <- lfc_mono25_withinComparision_ToTime0$log2FoldChange
lfc_mono_30withinComparisionToTime0.fc <- lfc_mono30_withinComparision_ToTime0$log2FoldChange
lfc_mono_35withinComparisionToTime0.fc <- lfc_mono35_withinComparision_ToTime0$log2FoldChange
lfc_mono_40withinComparisionToTime0.fc <- lfc_mono40_withinComparision_ToTime0$log2FoldChange
lfc_mono_45withinComparisionToTime0.fc <- lfc_mono45_withinComparision_ToTime0$log2FoldChange


save.image("RNAseq/LFC/outputFiles/PsyringaeLFC.RData")


##################These data were used for KEGG pathway analysis##################
###See "KeggPathwayAnalysis" script
