# KEGG pathway analysis results

The RNA-seq data consisted of a time course for each strain: An exponential phase time point (12.5h) followed by 5, stationary phase time points (25h, 30h, 35h, 40h, 45h).

For each strain, log2-fold change (LFC) values were calculated by comparing transcript counts at each stationary phase time point to transcript counts at the exponential phase time point.

These LFC changes were then mapped onto KEGG pathways using the pathview package in R.

To view KEGG pathways of interest, navigate into the strain directories to view the .png files. 

### Note, the default color scheme from pathview (Red: Increased expression, Green: Decreased expression) is opposite that for Figure 7 in the manuscript (Green: Increased expression, Red: Decreased expression)
