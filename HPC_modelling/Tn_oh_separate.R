library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)

setwd("~/Analysis/")

## Biomart ##
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")

#All_filters <- listAttributes(ensembl)
#View(All_filters)
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)


## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)

## prep data ##
print("preping data")
TPM_0h <- read.delim("~/Analysis/data/RNA/Tnaive_0h/quant.sf")
dim(TPM_0h)


dim(TPM_0h[TPM_0h$Name %in% tx2gene$ensembl_transcript_id_version,])

TPM_0h$Name <- mapply(strsplit(as.character(TPM_0h$Name),"\\."),FUN=function(x){(as.character(x)[1])})

dim(TPM_0h[TPM_0h$Name %in% Sequence_parameters_RNA$tx.id,])

TPM_0h <- data.frame("ID"=TPM_0h$Name,"TPM"=TPM_0h$TPM)
TPM_0h <- subset(TPM_0h, TPM_0h$TPM>0)
dim(TPM_0h)


#TPM_0h <- merge(TPM_0h, tx2gene, by.x="ID", by.y="ensembl_transcript_id", all.x=F)
dim(TPM_0h)
dim(TPM_0h[TPM_0h$transcript_biotype=="protein_coding",])
dim(TPM_0h[TPM_0h$ID %in% Sequence_parameters_RNA$tx.id,])
#dim(Sequence_parameters_RNA)
TPM_0h <- data.frame("ID"=TPM_0h$ID,"TPM"=TPM_0h$TPM)
TPM_0h <- merge(TPM_0h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
dim(TPM_0h)

print("integrating feature counts per gene")
registerDoMC(10)
TPM_0h_param <- ddply(TPM_0h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_0h_param$TPM <- NULL

print("integrating TPM per gene")
TPM_0h <- data.frame("ID"=TPM_0h$ensembl_gene_id ,"TPM"=TPM_0h$TPM)
TPM_0h <- ddply(TPM_0h,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_0h <- merge(TPM_0h,TPM_0h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_0h$TPM <- log10(TPM_0h$TPM)

TPM_0h[is.na(TPM_0h)]= 0
TPM_0h <- subset(TPM_0h, TPM_0h$TPM>=-1)
print(dim(TPM_0h)) # 13729  2682

rownames(TPM_0h) <- TPM_0h$ID
TPM_0h$ID <- NULL
TPM_0h$ensembl_gene_id <- NULL
TPM_0h$gene_biotype <- NULL
TPM_0h$gene_name <- NULL


write.table(TPM_0h,"TPM_0h_with_lib_per_gene.csv",row.names = T)


## Feature selection ##

## Tn 0h ##

print("removing useless columns")
TPM_0h[is.na(TPM_0h)]=0
TPM_0h <- TPM_0h[, -caret::nearZeroVar(TPM_0h, allowParallel = TRUE, uniqueCut = dim(TPM_0h)[1]*0.01)]
print(dim(TPM_0h)) # 13729  2637



gc()






## Training the RF model for Tn 0h ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_Tn_0h <- train(TPM~.,
                      data=TPM_0h,
                      method="rf",
                      ntree=1,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_0h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours


print(rf_TPM_Tn_0h)

print(rf_TPM_Tn_0h$finalModel)


saveRDS(rf_TPM_Tn_0h, "~/Analysis/models/rf_deep_5000trees_TPM_Tn_0h_29_06_2021.RDS")
