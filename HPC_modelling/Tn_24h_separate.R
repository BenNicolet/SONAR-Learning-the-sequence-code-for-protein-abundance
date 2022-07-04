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
TPM_24h <- read.delim("~/Analysis/data/RNA/Tnaive_act24h/quant.sf")
dim(TPM_24h)

TPM_24h <- data.frame("ID"=TPM_24h$Name,"TPM"=TPM_24h$TPM)
TPM_24h <- subset(TPM_24h, TPM_24h$TPM>0)
TPM_24h <- merge(TPM_24h, tx2gene, by.x="ID", by.y="ensembl_transcript_id_version", all.x=T)

registerDoMC(4)

TPM_24h <- ddply(TPM_24h,"external_gene_name", numcolwise(sum), .parallel = T, .progress = T)
TPM_24h$TPM <- log10(TPM_24h$TPM)
colnames(TPM_24h)[1]="ID"


TPM_24h <- merge(TPM_24h,Sequence_parameters_RNA,by="ID", all.x=T)
TPM_24h[is.na(TPM_24h)]= 0
TPM_24h <- subset(TPM_24h, TPM_24h$TPM>=-1)
print(dim(TPM_24h)) # 13729  2682

rownames(TPM_24h) <- TPM_24h$ID
TPM_24h$ID <- NULL


## Feature selection ##

registerDoMC(4)

## Tn 0h ##
TPM_24h[is.na(TPM_24h)]=0
TPM_24h <- TPM_24h[, -nearZeroVar(TPM_24h, allowParallel = T, uniqueCut = dim(TPM_24h)[1]*0.01)]
print(dim(TPM_24h)) # 13729  2682


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
rf_TPM_Tn_24h <- train(TPM~.,
                      data=TPM_24h,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_24h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours


print(rf_TPM_Tn_24h)

print(rf_TPM_Tn_24h$finalModel)


saveRDS(rf_TPM_Tn_24h, "~/Analysis/models/rf_TPM_Tn_24h_18_05_2021.RDS")
# rf_TPM_Tn_24h <- readRDS("~/Analysis/models/rf_TPM_Tn_24h_30_04_2021.RDS")
# 
# 
# imp_rf_TPM_Tn_24h<- data.frame(varImp(rf_TPM_Tn_24h$finalModel))
# imp_rf_TPM_Tn_24h$ID <- rownames(imp_rf_TPM_Tn_24h)

