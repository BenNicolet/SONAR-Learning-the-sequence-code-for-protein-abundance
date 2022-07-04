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
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Full_library_v1_with_protein_sites_30-03-21.csv", sep=" ", dec=".")


## prep data ##
CN_CD8_Tcm <- read.delim("~/Analysis/data/protein/CD8_Tcm_CN_log10.csv",sep=";",dec=",")
CN_CD8_Tcm$CN <- gsub(",",".",CN_CD8_Tcm$CN)
CN_CD8_Tcm$CN <- as.numeric(CN_CD8_Tcm$CN)
CN_CD8_Tcm <- subset(CN_CD8_Tcm,CN_CD8_Tcm$CN>0)
dim(CN_CD8_Tcm)

CN_CD8_Tcm <- merge(CN_CD8_Tcm,Sequence_parameters_protein,by="ID", all.x=T)
CN_CD8_Tcm[is.na(CN_CD8_Tcm)]= 0
print(dim(CN_CD8_Tcm)) # 8249 6803


rownames(CN_CD8_Tcm) <- CN_CD8_Tcm$ID
CN_CD8_Tcm$ID <- NULL


## Feature selection ##
registerDoMC(4)

CN_CD8_Tcm[is.na(CN_CD8_Tcm)]=0
CN_CD8_Tcm <- CN_CD8_Tcm[, -nearZeroVar(CN_CD8_Tcm, allowParallel = T, uniqueCut = dim(CN_CD8_Tcm)[1]*0.01)]
print(dim(CN_CD8_Tcm)) # 8249 2674


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
rf_CN_CD8_Tcm <- train(CN~.,
                      data=CN_CD8_Tcm,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Tcm)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_CD8_Tcm)

print(rf_CN_CD8_Tcm$finalModel)


saveRDS(rf_CN_CD8_Tcm, "~/Analysis/models/protein/rf_CN_CD8_Tcm_31_08_2021.RDS")




