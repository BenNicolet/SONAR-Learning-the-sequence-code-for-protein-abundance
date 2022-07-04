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
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")


## prep data ##
CN_CD8_Tn <- read.delim("~/Analysis/data/protein/CD8_Tn_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tn <- subset(CN_CD8_Tn,CN_CD8_Tn$CN>0)
dim(CN_CD8_Tn)

CN_CD8_Tn <- merge(CN_CD8_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tn[is.na(CN_CD8_Tn)]= 0
print(dim(CN_CD8_Tn)) # 7984 7126


rownames(CN_CD8_Tn) <- CN_CD8_Tn$ID
CN_CD8_Tn$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tn[is.na(CN_CD8_Tn)]=0
CN_CD8_Tn <- CN_CD8_Tn[, -nearZeroVar(CN_CD8_Tn, allowParallel = T, uniqueCut = dim(CN_CD8_Tn)[1]*0.01)]
print(dim(CN_CD8_Tn)) # 7984 2759

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
rf_CN_CD8_Tn <- train(CN~.,
                      data=CN_CD8_Tn,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Tn)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD8_Tn)
print(rf_CN_CD8_Tn$finalModel)

saveRDS(rf_CN_CD8_Tn, "~/Analysis/models/protein/rf_CN_CD8_Tn_09_09_2021.RDS")




