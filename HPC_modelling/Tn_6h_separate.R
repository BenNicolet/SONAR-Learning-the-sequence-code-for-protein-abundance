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

tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","external_gene_name","gene_biotype"), mart = ensembl)
tx2gene <- subset(tx2gene,tx2gene$gene_biotype=="protein_coding")


## parameters ##
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/Sequence_parameters_RNA_23-04-21.csv", sep=";", dec=",")


## prep data ##
TPM_6h <- read.delim("~/Analysis/data/RNA/Tnaive_act6h/quant.sf")
dim(TPM_6h)

TPM_6h <- data.frame("ID"=TPM_6h$Name,"TPM"=TPM_6h$TPM)
TPM_6h <- subset(TPM_6h, TPM_6h$TPM>0)
TPM_6h <- merge(TPM_6h, tx2gene, by.x="ID", by.y="ensembl_transcript_id_version", all.x=T)

registerDoMC(4)

TPM_6h <- ddply(TPM_6h,"external_gene_name", numcolwise(sum), .parallel = T, .progress = T)
TPM_6h$TPM <- log10(TPM_6h$TPM)
colnames(TPM_6h)[1]="ID"


TPM_6h <- merge(TPM_6h,Sequence_parameters_RNA,by="ID", all.x=T)
TPM_6h[is.na(TPM_6h)]= 0
TPM_6h <- subset(TPM_6h, TPM_6h$TPM>=-1)
print(dim(TPM_6h)) # 13729  2682

rownames(TPM_6h) <- TPM_6h$ID
TPM_6h$ID <- NULL


## Feature selection ##

registerDoMC(4)

## Tn 0h ##
TPM_6h[is.na(TPM_6h)]=0
TPM_6h <- TPM_6h[, -nearZeroVar(TPM_6h, allowParallel = T, uniqueCut = dim(TPM_6h)[1]*0.01)]
print(dim(TPM_6h)) # 13729  2682


gc()


## Training the RF model for Tn 0h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_Tn_6h <- train(TPM~.,
                      data=TPM_6h,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_6h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours


print(rf_TPM_Tn_6h)

print(rf_TPM_Tn_6h$finalModel)


saveRDS(rf_TPM_Tn_6h, "~/Analysis/models/rf_TPM_Tn_6h_18_05_2021.RDS")
# rf_TPM_Tn_6h <- readRDS("~/Analysis/models/rf_TPM_Tn_6h_30_04_2021.RDS")
# 
# 
# imp_rf_TPM_Tn_6h<- data.frame(varImp(rf_TPM_Tn_6h$finalModel))
# imp_rf_TPM_Tn_6h$ID <- rownames(imp_rf_TPM_Tn_6h)

