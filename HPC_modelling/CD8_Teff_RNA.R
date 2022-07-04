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
TPM_CD8_Teff <- read.delim("~/Analysis/data/RNA/CD8_Teff_TPM_log10.csv",sep=";",dec=",")
dim(TPM_CD8_Teff)

TPM_CD8_Teff <- merge(TPM_CD8_Teff,Sequence_parameters_RNA,by="ID", all.x=T)
TPM_CD8_Teff[is.na(TPM_CD8_Teff)]= 0
print(dim(TPM_CD8_Teff)) # 15225  6788

rownames(TPM_CD8_Teff) <- TPM_CD8_Teff$ID
TPM_CD8_Teff$ID <- NULL


## Feature selection ##
registerDoMC(4)

TPM_CD8_Teff[is.na(TPM_CD8_Teff)]=0
TPM_CD8_Teff <- TPM_CD8_Teff[, -nearZeroVar(TPM_CD8_Teff, allowParallel = T, uniqueCut = dim(TPM_CD8_Teff)[1]*0.01)]
print(dim(TPM_CD8_Teff)) # 15225  2635


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
rf_TPM_CD8_Teff <- train(TPM~.,
                      data=TPM_CD8_Teff,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_CD8_Teff)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_TPM_CD8_Teff)

print(rf_TPM_CD8_Teff$finalModel)


saveRDS(rf_TPM_CD8_Teff, "~/Analysis/models/rf_TPM_CD8_Teff_21_05_2021.RDS")

## 0h ##
# rf_TPM_CD8_Teff <- readRDS("~/Analysis/models/rf_TPM_CD8_Teff_18_05_2021.RDS")
# 
# imp_rf_TPM_CD8_Teff<- data.frame(varImp(rf_TPM_CD8_Teff$finalModel))
# imp_rf_TPM_CD8_Teff$ID <- rownames(imp_rf_TPM_CD8_Teff)

