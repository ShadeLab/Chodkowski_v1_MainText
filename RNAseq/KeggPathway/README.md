# KEGG pathway analysis results

The RNA-seq data consisted of a time course for each strain: An exponential phase time point (12.5h) followed by 5, stationary phase time points (25h, 30h, 35h, 40h, 45h).

For each strain, log2-fold change (LFC) values were calculated by comparing transcript counts at each stationary phase time point to transcript counts at the exponential phase time point.

These LFC changes were then mapped onto KEGG pathways using the pathview package in R.

To view KEGG pathways of interest, navigate back to the KeggPathway directory and then into the strain directories to view the .png files. 