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
TPM_CD4_Tem <- read.delim("~/Analysis/data/RNA/CD4_Tem_TPM_log10.csv",sep=";",dec=",")
dim(TPM_CD4_Tem)

TPM_CD4_Tem <- merge(TPM_CD4_Tem,Sequence_parameters_RNA,by="ID", all.x=T)
TPM_CD4_Tem[is.na(TPM_CD4_Tem)]= 0
print(dim(TPM_CD4_Tem)) # 15225  6788

rownames(TPM_CD4_Tem) <- TPM_CD4_Tem$ID
TPM_CD4_Tem$ID <- NULL


## Feature selection ##
registerDoMC(4)

TPM_CD4_Tem[is.na(TPM_CD4_Tem)]=0
TPM_CD4_Tem <- TPM_CD4_Tem[, -nearZeroVar(TPM_CD4_Tem, allowParallel = T, uniqueCut = dim(TPM_CD4_Tem)[1]*0.01)]
print(dim(TPM_CD4_Tem)) # 15225  2635


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
rf_TPM_CD4_Tem <- train(TPM~.,
                      data=TPM_CD4_Tem,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_CD4_Tem)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_TPM_CD4_Tem)

print(rf_TPM_CD4_Tem$finalModel)


saveRDS(rf_TPM_CD4_Tem, "~/Analysis/models/rf_TPM_CD4_Tem_21_05_2021.RDS")

## 0h ##
# rf_TPM_CD4_Tem <- readRDS("~/Analysis/models/rf_TPM_CD4_Tem_18_05_2021.RDS")
# 
# imp_rf_TPM_CD4_Tem<- data.frame(varImp(rf_TPM_CD4_Tem$finalModel))
# imp_rf_TPM_CD4_Tem$ID <- rownames(imp_rf_TPM_CD4_Tem)

