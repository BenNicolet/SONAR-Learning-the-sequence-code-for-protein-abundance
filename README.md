# Learning the sequence code for mRNA and protein abundance with SONAR 

We developed SONAR an auditable machine learning pipeline to understand the relation between sequence features (e.g. RNA-binding protein motifs, GC-content,...) and mRNA and protein abundance in immune cells and cell lines. 
Our preprint is available here : https://www.biorxiv.org/content/10.1101/2023.09.01.555843v1

Below, we describe the steps needed to reproduce our findings, or adapt it to your data. 



### 1. Preparing the sequence feature library
Here we prepare the sequence feature library. We used many feature types (see Fig S1 for all details), including RNA motifs, GC content, codon and amino-acid usage. You can adapt these as you wish (just make sure these are either categorical or continuous). 
First, we acquire the 5'UTR, CDS, and 3'UTR fasta sequences from Ensembl Biomart (https://www.ensembl.org/biomart/martview), and we feed these to the R-script     “./1_Library_preparation/RNA&protein_library_V2_generation.R” . 
The RBP motifs were obtained from https://attract.cnic.es/index and are referred to in scripts as "ATtRACT_db.txt”. 
The library generation scripts will search for sequence features in the fasta sequences in the entire (coding) transcriptome. You can find more details in the script’s annotation. Of note, we generated 2 libraries 1 for mRNA and 1 for protein features and prediction. 


### 2. Getting mRNA and protein abundance data 
Here we used previously published data (see Table S1 of our preprint). The data are available in supplementary tables and in the "input data" folder on GitHub. We recommend using TPM (Transcript per kilobase million; computed with Salmon) as it allows a fairer inter-sample comparison. For Mass Spectrometry-based proteomics data, we recommend using either iBAQ or protein copy number (CN; we computed CN from LFQ using the proteomic ruler package of Perseus). 
Once the data are prepared, we log10-transform and aggregate these in 1 table to facilitate the modelling in batch.


### 3. Modelling the mRNA and protein abundance from sequence features
Here we model the abundance data using the sequence feature library. To facilitate the modelling of many samples, we loop the modelling part of the script (see “3. Modelling” folder and scripts for more details). 	
In brief, the scripts will A) get the data and merge these to the sequence feature library, B) define the modelling parameters, C) split the data (at random) in a test and training set (preventing any data leakage) and save these, D) model (using XGBoost) and save the model, E) use the held-out test set and the model to evaluate the prediction accuracy, F) extract the feature importance (aka measure of feature contribution to the model) and save these. 


### 4. Data visualization 
This part is not directly part of SONAR, and it is “work in progress” and should be used as inspiration for data analysis. 
