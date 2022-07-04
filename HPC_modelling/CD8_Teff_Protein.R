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
CN_CD8_Teff <- read.delim("~/Analysis/data/protein/CD8_Teff_CN_log10.csv",sep=";",dec=",")
CN_CD8_Teff$CN <- gsub(",",".",CN_CD8_Teff$CN)
CN_CD8_Teff$CN <- as.numeric(CN_CD8_Teff$CN)
CN_CD8_Teff <- subset(CN_CD8_Teff,CN_CD8_Teff$CN>0)
dim(CN_CD8_Teff)

CN_CD8_Teff <- merge(CN_CD8_Teff,Sequence_parameters_protein,by="ID", all.x=T)
CN_CD8_Teff[is.na(CN_CD8_Teff)]= 0
print(dim(CN_CD8_Teff)) # 8190 6803


rownames(CN_CD8_Teff) <- CN_CD8_Teff$ID
CN_CD8_Teff$ID <- NULL


## Feature selection ##
registerDoMC(4)

CN_CD8_Teff[is.na(CN_CD8_Teff)]=0
CN_CD8_Teff <- CN_CD8_Teff[, -nearZeroVar(CN_CD8_Teff, allowParallel = T, uniqueCut = dim(CN_CD8_Teff)[1]*0.01)]
print(dim(CN_CD8_Teff)) # 8190 2680


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
rf_CN_CD8_Teff <- train(CN~.,
                      data=CN_CD8_Teff,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Teff)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_CD8_Teff)

print(rf_CN_CD8_Teff$finalModel)


saveRDS(rf_CN_CD8_Teff, "~/Analysis/models/protein/rf_CN_CD8_Teff_31_08_2021.RDS")




