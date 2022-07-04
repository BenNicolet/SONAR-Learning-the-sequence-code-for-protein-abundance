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


## parameters ##
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/Sequence_parameters_RNA_23-04-21.csv", sep=";", dec=",")


## prep data ##
TPM_CD8_Tn <- read.delim("~/Analysis/data/RNA/CD8_Tn_TPM_log10.csv",sep=";",dec=",")
dim(TPM_CD8_Tn)

TPM_CD8_Tn <- merge(TPM_CD8_Tn,Sequence_parameters_RNA,by="ID", all.x=T)
TPM_CD8_Tn[is.na(TPM_CD8_Tn)]= 0
print(dim(TPM_CD8_Tn)) # 15225  6788

rownames(TPM_CD8_Tn) <- TPM_CD8_Tn$ID
TPM_CD8_Tn$ID <- NULL


## Feature selection ##
registerDoMC(4)

TPM_CD8_Tn[is.na(TPM_CD8_Tn)]=0
TPM_CD8_Tn <- TPM_CD8_Tn[, -nearZeroVar(TPM_CD8_Tn, allowParallel = T, uniqueCut = dim(TPM_CD8_Tn)[1]*0.01)]
print(dim(TPM_CD8_Tn)) # 15225  2635


gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_CD8_Tn <- train(TPM~.,
                      data=TPM_CD8_Tn,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_CD8_Tn)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_TPM_CD8_Tn)

print(rf_TPM_CD8_Tn$finalModel)


saveRDS(rf_TPM_CD8_Tn, "~/Analysis/models/rf_TPM_CD8_Tn_21_05_2021.RDS")

## 0h ##
# rf_TPM_CD8_Tn <- readRDS("~/Analysis/models/rf_TPM_CD8_Tn_18_05_2021.RDS")
# 
# imp_rf_TPM_CD8_Tn<- data.frame(varImp(rf_TPM_CD8_Tn$finalModel))
# imp_rf_TPM_CD8_Tn$ID <- rownames(imp_rf_TPM_CD8_Tn)

