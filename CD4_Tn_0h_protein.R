library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)

setwd("/home/ben/Analysis/RF_human/T_cell_subsets/")

## parameters ##
print("importing library")
#Sequence_parameters_Geiger <- read.delim("/home/ben/Analysis/RF_human/sequences/Sequence_parameters_for_Geiger_data_prot_libv1_21_05_21.csv", sep=";", dec=",")
Sequence_parameters_Geiger <- read.delim("/home/ben/Analysis/RF_human/Library_V2_Aug2021/Protein_per_gene_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv",sep = "\t")

Sequence_parameters_Geiger$tx.id <- NULL
Sequence_parameters_Geiger$Entry <- NULL
Sequence_parameters_Geiger$ensembl_gene_id.x <- NULL


## prep data ##
print("importing and preparing data")
CN_0h <- read.delim("/home/ben/Analysis/RF_human/T_cell_subsets/Protein/CD4_Tn_0h_CN_log10.csv",sep = ";",dec = ",")
CN_0h <- merge(CN_0h,Sequence_parameters_Geiger,by.x="ID",by.y="gene_name",all.x=F)
print(dim(CN_0h)) # 6993 6803


rownames(CN_0h) <- CN_0h$ID
CN_0h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_0h[is.na(CN_0h)]=0
CN_0h <- CN_0h[, -nearZeroVar(CN_0h, allowParallel = T, uniqueCut = dim(CN_0h)[1]*0.01)]
print(dim(CN_0h)) # 13729  2682

gc()


## Training the RF model for Tn 0h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h_ns5 <- train(CN~.,
                      data=CN_0h,
                      method="rf",
                      ntree=1000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_Tn_0h_ns5)

print(rf_CN_Tn_0h_ns5$finalModel)


saveRDS(rf_CN_Tn_0h_ns5, "/home/ben/Analysis/RF_human/WIP_models/rf_CN_Tn_0h_Lib_V2_nodesize_5_03092021.RDS")




## Training the RF model for Tn 0h nodesize 10 ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h_ns10 <- train(CN~.,
                     data=CN_0h,
                     method="rf",
                     ntree=1000,
                     nodesize=10,
                     trControl=control,
                     metric="Rsquared",
                     tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                     na.action = na.omit,
                     importance = TRUE,
                     verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_Tn_0h_ns10)

print(rf_CN_Tn_0h_ns10$finalModel)


saveRDS(rf_CN_Tn_0h_ns10, "/home/ben/Analysis/RF_human/WIP_models/rf_CN_Tn_0h_Lib_V2_nodesize_10_03092021.RDS")



## Training the RF model for Tn 0h nodesize 20 ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h_ns20 <- train(CN~.,
                          data=CN_0h,
                          method="rf",
                          ntree=1000,
                          nodesize=20,
                          trControl=control,
                          metric="Rsquared",
                          tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                          na.action = na.omit,
                          importance = TRUE,
                          verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_Tn_0h_ns20)

print(rf_CN_Tn_0h_ns20$finalModel)


saveRDS(rf_CN_Tn_0h_ns20, "/home/ben/Analysis/RF_human/WIP_models/rf_CN_Tn_0h_Lib_V2_nodesize_20_03092021.RDS")




## Training the RF model for Tn 0h nodesize 30 ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h_ns30 <- train(CN~.,
                          data=CN_0h,
                          method="rf",
                          ntree=1000,
                          nodesize=30,
                          trControl=control,
                          metric="Rsquared",
                          tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                          na.action = na.omit,
                          importance = TRUE,
                          verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_Tn_0h_ns30)

print(rf_CN_Tn_0h_ns30$finalModel)


saveRDS(rf_CN_Tn_0h_ns30, "/home/ben/Analysis/RF_human/WIP_models/rf_CN_Tn_0h_Lib_V2_nodesize_30_03092021.RDS")



## Training the RF model for Tn 0h nodesize 30 ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h_ns50 <- train(CN~.,
                          data=CN_0h,
                          method="rf",
                          ntree=1000,
                          nodesize=50,
                          trControl=control,
                          metric="Rsquared",
                          tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                          na.action = na.omit,
                          importance = TRUE,
                          verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 


print(rf_CN_Tn_0h_ns50)

print(rf_CN_Tn_0h_ns50$finalModel)


saveRDS(rf_CN_Tn_0h_ns50, "/home/ben/Analysis/RF_human/WIP_models/rf_CN_Tn_0h_Lib_V2_nodesize_50_03092021.RDS")


