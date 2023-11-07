# Learning the sequence code for mRNA and protein abundance with SONAR 

We develloped SONAR an auditable machine learning pipeline to understand the relation between sequence features (e.g. RNA-binding protein motifs, GC-content,...) and mRNA and protein abundance in immune cells and cell lines. 
Our preprint is available here : https://www.biorxiv.org/content/10.1101/2023.09.01.555843v1

Below, we describe the steps needed to reproduce our findings, or adapt it to your data. 

### Getting mRNA and protein abundance data 
Here we used previsouly published data (see Fig S1 of our preprint). The data are available in supplementary tables and in the "input data" folder on github. We recommend using TPM (Transcript per kilobase million; computed with Salmon) as it allows a fairer inter-sample comparison. For Mass Spectrometry-based proteomics data, we recommend using either iBAQ or protein copy number (CN; we computed CN from LFQ using the proteomic ruler package of Perseus). 

Once the data are prepared, we log10-transform and aggregate these in 1 table to facilitate the modelling in batch.

### preparing the sequence feature library
