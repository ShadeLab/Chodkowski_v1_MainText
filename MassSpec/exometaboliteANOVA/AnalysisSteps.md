# Step 1: Obtain and export only L1 exometabolites

This step can be found in the ReleasedExometabolites_PolarPos.R and ReleasedExometabolites_PolarNeg.R scripts

Specifically under the section: Identification of released metabolites
```
write.csv(allMSL1.O,"MassSpec/releaseAnalysis/MSMS/outputFiles/PolarPos_JGI_MSL1_IDs.csv",row.names=FALSE)
#Export label names
write.csv(allJGI,"MassSpec/releaseAnalysis/MSMS/outputFiles/PolarPos_JGI_MSL1_Labels.csv",row.names=FALSE)
```
This was used to generate files: <br />
1) PolarPos_L1_exometabolites.csv <br />
2) PolarNeg_L1_exometabolites.csv

# Step 2: Normalize dataset in MetaboAnalyst

Files 1 & 2 from Step 1 were uploaded to MetaboAnalyst

Datasets were normalized by a reference feature(13C-15N-Alanine for Polar Positive and 13C-15N-Proline for Polar Negative) and cubed-root trasnformed.

Normalized datasets were exported to generate files: <br />

1) PolarPos_L1_exometabolites_normalized.csv <br />
2) PolarNeg_L1_exometabolites_normalized.csv

# Step 3: Combine datasets from both ionization modes

Exometabolites in files 1 & 2 from step 2 were combined to generate: <br />
L1_exometabolites_combinedIonizations_normalized.csv

# Step 4: Upload to MetaboAnalyst and perform ANOVA.

The result was downloaded and generated file: <br />
ANOVA_results.csv

Specifically, this file was used to focus on ANOVA results of the most accumulated exometabolite in each strain. 
