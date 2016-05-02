# Cohens-d-replication
This repository contains all files that are needed to reproduce the findings that were presented in the paper "Replication of Cohen's d: Quantifying Evidence using the Bayes Factor". Please contact Millitza D. Kroonenberg (millitza@gmail.com) or prof. dr. Herbert J. A. Hoijtink (h.hoijtink@uu.nl) if you have any questions related to this repository. The repository contains the following files:
## function.R
This R file contains annotated code that can be used to replicate the findings in the paper. The code contains the two functions as described in manualRfunction.pdf. For users who are unfamiliar with R there is an application with user interface available (see folder Application).
## manualRfunction.pdf
This manual contains information on the functions that were used in the paper and can be found in function.R. Both the required input and returned objects are discussed.
## manuscript.pdf
This is the manuscript of the paper "Replication of Cohen's d: Quantifying Evidence using the Bayes Factor". This version was handed in as the thesis for the Utrecht University master programme Methodology and Statistics for the Behavioral, Biomedical, and Social Sciences.
## reprod.R
This R file can be run to acquire the exact results as presented in the paper (sections "Evidence for replication of Cohen's d" and "Data example: Mindset conditions and ambivalence"). The seed numbers that were used are also specified within this function file.

## Application (folder)
This folder contains the files that are needed to calculate Bayes factors that quantify effect size replication evidence in an online application with user interface. Note that you must save the folder "Application" and store the file "function.R" in the root-folder.
### manualShiny.pdf
This manual describes how to use the Shiny application in more detail. Users that are more advanced in using R or want more options are encouraged to use the R function.
### ui.R
This file contains the code that is needed to create a user interface in Shiny.
### server.R
This file contains the code that is needed to calculate the components of the Shiny application.

## BIEMS data (folder)
This folder contains the data files as generated by BIEMS, using seed=864 and "exact" is enabled. Each row represents a case/participant, while the first column contains its value on the variable of interest and the second column indicates group membership (1 or 2). The file name is created as follows: NXX represents the number of cases/participants per group (N20 indicates 20 cases/participants per group), dXX represents the group mean difference while the Group 1 mean is zero per definition (d02 indicates Group 2 mean is 0.2).
### BIEMSdataPrep.R
This R file contains a short script that can be used to load BIEMS data and transform to a data file that can be read by the functions in function.R.
### BIEMS program (sub folder)
- data.txt
The resulting data file as generated by BIEMS.
- genmvldata.exe
The program that BIEMS uses to generate (exact) data sets. "genmvldata" is part of the official BIEMS program (downloaded from http://informative-hypotheses.sites.uu.nl/software/biems/). This program was used to generate the data sets described in the section "Evidence for replication of Cohen's d" and the original data set dat was used in the section "Data example".
- input.txt
This is the input file that defines the properties of the generated data set. This file is generated automatically by BIEMS.
- manual-genmvldata.pdf
User manual for the data generating program genmvldata.exe.

## Data example (folder)
### Original study (sub folder)
##### data_manipulation_orig.R
This R-file contains all the code that is required to transform the three data sets that were generated in BIEMS (data_neutral.txt, data_onesided.txt, data_twosided.txt) into a format that can be analyzed by the R-functions in function.R and reprod.R.
##### data_neutral.txt, data_onesided.txt, data_twosided.txt
These data sets were generated using BIEMS (seed=4578). The data is based on the summary statistics as provided by Henderson et al. (2008) on the original data set. The data sets match the summary statistics exactly since the option to generate exact data sets was enabled in BIEMS. The files data_neutral.txt, data_onesided.txt, data_twosided.txt contain generated data sets for the neutral, one-sided, and two-sided mindset conditions of the original study respectively.
#### Data_orig_COND1_COND2.txt
These are the three data files that can be analyzed directly using reprod.R. These files were created using data manipulation orig.R.
Each row is a participant. Explanation of variables in the data set:
[,1]	score on outcome variable (ambivalence)
[,2]	membership of COND1 (0=no, 1=yes)
[,3]	membership of COND2 (0=no, 1=yes)
#### Henderson_et_al_2008.pdf
This is the published journal article that was used as an example of original data in the manuscript (section "Data example").
### Replication study (sub folder)
#### data_manipulation_repl.R
This R-file contains all the code that is required to transform the data set that was provided by the researchers of the replication study (reconciling.for.repro.txt) into a format that can be analyzed by the R-functions in function.R and reprod.R. The first data manipulation steps were suggested by Lane and Gazarian in order to match the data editing process in the original study.
#### Data_repl_COND1_COND2.txt
These are the three data files that can be analyzed directly using reprod.R. These files were created using data manipulation orig.R.
Each row is a participant. Explanation of variables in the data set:
[,1]	score on outcome variable (ambivalence)
[,2]	membership of COND1 (0=no, 1=yes)
[,3]	membership of COND2 (0=no, 1=yes)
#### henderson.replication.report.final.pdf
This is the official research report accompanying the replication attempt and data by Lane and Gazarian (2015).
#### reconciling.for.repro.txt
This data file was downloaded via https://osf.io/3xv8u/ and contained the data set of the replication study.

## Figures (folder)
This folder contains all files that were used to create the figures as presented in the paper.
