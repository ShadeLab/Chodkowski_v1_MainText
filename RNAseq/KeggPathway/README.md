# KEGG pathway analysis results

The RNA-seq data consisted of a time course for each strain: An exponential phase time point (12.5) followed by 5, stationary phase time points (25, 30, 35, 40, 45h).

For all strains, log2-fold change (LFC) values were calculated for transcript counts at each stationary phase time point by comparing transcript counts to the exponential phase time point.

These LFC changes were then mapped onto KEGG pathways using the pathview package in R.

To view mapped pathways of interest, navigate into the strain directories and open the .png files. 
