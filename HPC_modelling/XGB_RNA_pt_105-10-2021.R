library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)
library(xgboost)
library(ggpointdensity)
library(viridis)




###__________________________________________________________________________________________###
###____________________________________________Ex vivo_______________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")

## Biomart ##
# ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
# tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)
#
#
# ## parameters ##
# print("importing lib")
# Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
# dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###____________________________________________Tn_0h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for Tn 0h")
# TPM_0h <- read.delim("~/Analysis/data/RNA/Tnaive_0h/quant.sf")
#
# TPM_0h$Name <- mapply(strsplit(as.character(TPM_0h$Name),"\\."),FUN=function(x){(as.character(x)[1])})
#
# TPM_0h <- data.frame("ID"=TPM_0h$Name,"TPM"=TPM_0h$TPM)
# TPM_0h <- subset(TPM_0h, TPM_0h$TPM>0)
# TPM_0h <- data.frame("ID"=TPM_0h$ID,"TPM"=TPM_0h$TPM)
# TPM_0h <- merge(TPM_0h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
#
#
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_0h_param <- ddply(TPM_0h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_0h_param$TPM <- NULL
#
# print("integrating TPM per gene")
# TPM_0h <- data.frame("ID"=TPM_0h$ensembl_gene_id ,"TPM"=TPM_0h$TPM)
# TPM_0h <- ddply(TPM_0h,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_0h <- merge(TPM_0h,TPM_0h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
#
#
# TPM_0h$TPM <- log10(TPM_0h$TPM)
#
# TPM_0h[is.na(TPM_0h)]= 0
# TPM_0h <- subset(TPM_0h, TPM_0h$TPM>=-1)
# print(dim(TPM_0h)) # 13729  2682
#
# rownames(TPM_0h) <- TPM_0h$ID
# TPM_0h$ID <- NULL
# TPM_0h$ensembl_gene_id <- NULL
# TPM_0h$gene_biotype <- NULL
# TPM_0h$gene_name <- NULL
#
#
# #write.table(TPM_0h,"~/Analysis/counts&libs/TPM_0h_with_lib_per_gene.csv",row.names = T)
#
#
# ## Feature selection ##
# ## Tn 0h ##
#
# print("removing useless columns")
# TPM_0h[is.na(TPM_0h)]=0
# TPM_0h <- TPM_0h[, -caret::nearZeroVar(TPM_0h, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_0h)) # 13729  2637

# gc()
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_0h), 0.8 * nrow(TPM_0h))
# TPM_0h_train <- TPM_0h[train_row,]
# TPM_0h_test <- TPM_0h[-train_row,]
#
# dim(TPM_0h_train)
# dim(TPM_0h_test)
#
# write.table(TPM_0h_train,"~/Analysis/split_dataset/TPM_0h_train",sep=";")
# write.table(TPM_0h_test,"~/Analysis/split_dataset/TPM_0h_test",sep=";")
#
# TPM_0h_train <- read.delim("~/Analysis/split_dataset/TPM_0h_train",sep=";")
# TPM_0h_test <- read.delim("~/Analysis/split_dataset/TPM_0h_test",sep=";")
#
#
#
# ## Training the xgb model for Tn 0h ##
# print("modeling")
# ## Training the XGB model for Tn 0h ##
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_Tn_0h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                              max_depth = 6,
#                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                              ## The values below are default values in the sklearn-api.
#                              eta = 0.05,
#                              gamma=0,
#                              min_child_weight = 0.9,
#                              subsample = 1)
#
# start_time <- Sys.time()
# xgb_TPM_Tn_0h <- train(TPM~.,
#                         data=TPM_0h_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_Tn_0h,
#                         na.action = na.omit,
#                         nthread=72,
#                         verbose = TRUE)
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of 16.78211 hours
#
# print(xgb_TPM_Tn_0h)
#
#
# saveRDS(xgb_TPM_Tn_0h, "~/Analysis/models/xgb_1000nrounds_TPM_Tn_0h_15_10_2021.RDS")
#
# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_1000nrounds_TPM_Tn_0h_15_10_2021.RDS")
#
xgb_TPM_Tn_0h_predict_test <- predict(xgb_TPM_Tn_0h,TPM_0h_test)
xgb_TPM_Tn_0h_predict_test <- data.frame(xgb_TPM_Tn_0h_predict_test)
print(lm(xgb_TPM_Tn_0h_predict_test$xgb_TPM_Tn_0h_predict_test~TPM_0h_test$TPM))
print(cor(xgb_TPM_Tn_0h_predict_test$xgb_TPM_Tn_0h_predict_test,TPM_0h_test$TPM, method = "pearson", use = "complete.obs")^2)
#
#
# gc()
# rm(TPM_0h,TPM_0h_param)






###__________________________________________________________________________________________###
###____________________________________________Tn_6h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data")
# TPM_6h <- read.delim("~/Analysis/data/RNA/Tnaive_act6h/quant.sf")
#
# TPM_6h$Name <- mapply(strsplit(as.character(TPM_6h$Name),"\\."),FUN=function(x){(as.character(x)[1])})
#
# TPM_6h <- data.frame("ID"=TPM_6h$Name,"TPM"=TPM_6h$TPM)
# TPM_6h <- subset(TPM_6h, TPM_6h$TPM>0)
# TPM_6h <- data.frame("ID"=TPM_6h$ID,"TPM"=TPM_6h$TPM)
# TPM_6h <- merge(TPM_6h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
#
#
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_6h_param <- ddply(TPM_6h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_6h_param$TPM <- NULL
#
# print("integrating TPM per gene")
# TPM_6h <- data.frame("ID"=TPM_6h$ensembl_gene_id ,"TPM"=TPM_6h$TPM)
# TPM_6h <- ddply(TPM_6h,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_6h <- merge(TPM_6h,TPM_6h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
#
#
# TPM_6h$TPM <- log10(TPM_6h$TPM)
#
# TPM_6h[is.na(TPM_6h)]= 0
# TPM_6h <- subset(TPM_6h, TPM_6h$TPM>=-1)
# print(dim(TPM_6h)) # 13729  2682
#
# rownames(TPM_6h) <- TPM_6h$ID
# TPM_6h$ID <- NULL
# TPM_6h$ensembl_gene_id <- NULL
# TPM_6h$gene_biotype <- NULL
# TPM_6h$gene_name <- NULL
#
#
# #write.table(TPM_6h,"~/Analysis/counts&libs/TPM_6h_with_lib_per_gene.csv",row.names = T)
#
#
# ## Feature selection ##
# ## Tn 6h ##
#
# print("removing useless columns")
# TPM_6h[is.na(TPM_6h)]=0
# TPM_6h <- TPM_6h[, -caret::nearZeroVar(TPM_6h, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_6h)) # 13729  2637
#
# gc()
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_6h), 0.8 * nrow(TPM_6h))
# TPM_6h_train <- TPM_6h[train_row,]
# TPM_6h_test <- TPM_6h[-train_row,]
#
# dim(TPM_6h_train)
# dim(TPM_6h_test)
#
# write.table(TPM_6h_train,"~/Analysis/split_dataset/TPM_6h_train",sep=";")
# write.table(TPM_6h_test,"~/Analysis/split_dataset/TPM_6h_test",sep=";")

TPM_6h_train <- read.delim("~/Analysis/split_dataset/TPM_6h_train",sep=";")
TPM_6h_test <- read.delim("~/Analysis/split_dataset/TPM_6h_test",sep=";")
#
#
#
# ## Training the xgb model for Tn 6h ##
# print("modeling")
# ## Training the XGB model for Tn 6h ##
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_Tn_6h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                              max_depth = 6,
#                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                              ## The values below are default values in the sklearn-api.
#                              eta = 0.05,
#                              gamma=0,
#                              min_child_weight = 0.9,
#                              subsample = 1)
#
# start_time <- Sys.time()
# xgb_TPM_Tn_6h <- train(TPM~.,
#                        data=TPM_6h_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_Tn_6h,
#                        na.action = na.omit,
#                        nthread=72,
#                        verbose = TRUE)
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of 16.78211 hours
#
# print(xgb_TPM_Tn_6h)
#
# saveRDS(xgb_TPM_Tn_6h, "~/Analysis/models/xgb_1000nrounds_TPM_Tn_6h_15_10_2021.RDS")
xgb_TPM_Tn_6h <- readRDS("~/Analysis/models/xgb_1000nrounds_TPM_Tn_6h_15_10_2021.RDS")
#
xgb_TPM_Tn_6h_predict_test <- predict(xgb_TPM_Tn_6h,TPM_6h_test)
xgb_TPM_Tn_6h_predict_test <- data.frame(xgb_TPM_Tn_6h_predict_test)
print(lm(xgb_TPM_Tn_6h_predict_test$xgb_TPM_Tn_6h_predict_test~TPM_6h_test$TPM))
print(cor(xgb_TPM_Tn_6h_predict_test$xgb_TPM_Tn_6h_predict_test,TPM_6h_test$TPM, method = "pearson", use = "complete.obs")^2)

#
#
#
#
# gc()
#
# rm(TPM_6h,TPM_6h_param)







###__________________________________________________________________________________________###
###____________________________________________Tn_24h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data")
# TPM_24h <- read.delim("~/Analysis/data/RNA/Tnaive_act24h/quant.sf")
#
# TPM_24h$Name <- mapply(strsplit(as.character(TPM_24h$Name),"\\."),FUN=function(x){(as.character(x)[1])})
#
# TPM_24h <- data.frame("ID"=TPM_24h$Name,"TPM"=TPM_24h$TPM)
# TPM_24h <- subset(TPM_24h, TPM_24h$TPM>0)
# TPM_24h <- data.frame("ID"=TPM_24h$ID,"TPM"=TPM_24h$TPM)
# TPM_24h <- merge(TPM_24h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
#
#
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_24h_param <- ddply(TPM_24h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_24h_param$TPM <- NULL
#
# print("integrating TPM per gene")
# TPM_24h <- data.frame("ID"=TPM_24h$ensembl_gene_id ,"TPM"=TPM_24h$TPM)
# TPM_24h <- ddply(TPM_24h,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_24h <- merge(TPM_24h,TPM_24h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
#
#
# TPM_24h$TPM <- log10(TPM_24h$TPM)
#
# TPM_24h[is.na(TPM_24h)]= 0
# TPM_24h <- subset(TPM_24h, TPM_24h$TPM>=-1)
# print(dim(TPM_24h)) # 13729  2682
#
# rownames(TPM_24h) <- TPM_24h$ID
# TPM_24h$ID <- NULL
# TPM_24h$ensembl_gene_id <- NULL
# TPM_24h$gene_biotype <- NULL
# TPM_24h$gene_name <- NULL
#
#
# #write.table(TPM_24h,"~/Analysis/counts&libs/TPM_24h_with_lib_per_gene.csv",row.names = T)
#
#
# ## Feature selection ##
# ## Tn 24h ##
#
# print("removing useless columns")
# TPM_24h[is.na(TPM_24h)]=0
# TPM_24h <- TPM_24h[, -caret::nearZeroVar(TPM_24h, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_24h)) # 13729  2637
#
# gc()
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_24h), 0.8 * nrow(TPM_24h))
# TPM_24h_train <- TPM_24h[train_row,]
# TPM_24h_test <- TPM_24h[-train_row,]
#
# dim(TPM_24h_train)
# dim(TPM_24h_test)
#
# write.table(TPM_24h_train,"~/Analysis/split_dataset/TPM_24h_train",sep=";")
# write.table(TPM_24h_test,"~/Analysis/split_dataset/TPM_24h_test",sep=";")

TPM_24h_train <- read.delim("~/Analysis/split_dataset/TPM_24h_train",sep=";")
TPM_24h_test <- read.delim("~/Analysis/split_dataset/TPM_24h_test",sep=";")
#
#
#
# ## Training the xgb model for Tn 24h ##
# print("modeling")
# ## Training the XGB model for Tn 24h ##
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_Tn_24h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                              max_depth = 6,
#                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                              ## The values below are default values in the sklearn-api.
#                              eta = 0.05,
#                              gamma=0,
#                              min_child_weight = 0.9,
#                              subsample = 1)
#
# start_time <- Sys.time()
# xgb_TPM_Tn_24h <- train(TPM~.,
#                        data=TPM_24h_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_Tn_24h,
#                        na.action = na.omit,
#                        nthread=72,
#                        verbose = TRUE)
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of 16.78211 hours
#
# print(xgb_TPM_Tn_24h)
#
#
# saveRDS(xgb_TPM_Tn_24h, "~/Analysis/models/xgb_1000nrounds_TPM_Tn_24h_15_10_2021.RDS")
xgb_TPM_Tn_24h <- readRDS( "~/Analysis/models/xgb_1000nrounds_TPM_Tn_24h_15_10_2021.RDS")

xgb_TPM_Tn_24h_predict_test <- predict(xgb_TPM_Tn_24h,TPM_24h_test)
xgb_TPM_Tn_24h_predict_test <- data.frame(xgb_TPM_Tn_24h_predict_test)
print(lm(xgb_TPM_Tn_24h_predict_test$xgb_TPM_Tn_24h_predict_test~TPM_24h_test$TPM))
print(cor(xgb_TPM_Tn_24h_predict_test$xgb_TPM_Tn_24h_predict_test,TPM_24h_test$TPM, method = "pearson", use = "complete.obs")^2)
#
#
#
# gc()
# rm(TPM_24h,TPM_24h_param)
#









###__________________________________________________________________________________________###
###____________________________________________ex vivo________________________________________###
###__________________________________________________________________________________________###


## parameters ##
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")



###__________________________________________________________________________________________###
###____________________________________________CD8_Tn________________________________________###
###__________________________________________________________________________________________###


## prep data ##
# CN_CD8_Tn <- read.delim("~/Analysis/data/protein/CD8_Tn_CN_log10.csv",sep=";",dec = ",")
# CN_CD8_Tn <- subset(CN_CD8_Tn,CN_CD8_Tn$CN>0)
# dim(CN_CD8_Tn)
#
# CN_CD8_Tn <- merge(CN_CD8_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_CD8_Tn[is.na(CN_CD8_Tn)]= 0
# print(dim(CN_CD8_Tn)) # 7984 7126
#
#
# rownames(CN_CD8_Tn) <- CN_CD8_Tn$ID
# CN_CD8_Tn$ID <- NULL
#
#
# ## Feature selection ##
# registerDoMC(6)
#
# CN_CD8_Tn[is.na(CN_CD8_Tn)]=0
# CN_CD8_Tn <- CN_CD8_Tn[, -nearZeroVar(CN_CD8_Tn, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_CD8_Tn)) # 7984 2759
#
# gc()
#
#
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_CD8_Tn), 0.8 * nrow(CN_CD8_Tn))
# CN_CD8_Tn_train <- CN_CD8_Tn[train_row,]
# CN_CD8_Tn_test <- CN_CD8_Tn[-train_row,]
#
# dim(CN_CD8_Tn_train)
# dim(CN_CD8_Tn_test)
#
# write.table(CN_CD8_Tn_train,"~/Analysis/split_dataset/CN_CD8_Tn_train",sep=";")
# write.table(CN_CD8_Tn_test,"~/Analysis/split_dataset/CN_CD8_Tn_test",sep=";")

CN_CD8_Tn_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tn_train",sep=";")
CN_CD8_Tn_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tn_test",sep=";")
#
#
#
# ## Training the xgb model ##
# print("modeling CD8 TN CN")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_CD8_TN_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                               max_depth = 6,
#                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                               ## The values below are default values in the sklearn-api.
#                               eta = 0.05,
#                               gamma=0,
#                               min_child_weight = 0.9,
#                               subsample = 1)
#
# start_time <- Sys.time()
# xgb_CD8_TN_CN <- train(CN~.,
#                         data=CN_CD8_Tn_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_CD8_TN_CN,
#                         na.action = na.omit,
#                         nthread=72,
#                         verbose = TRUE)
#
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
#
# saveRDS(xgb_CD8_TN_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tn_18_10_2021.RDS")
xgb_CD8_TN_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tn_18_10_2021.RDS")
#
# print(xgb_CD8_TN_CN)

xgb_CD8_TN_CN_predict_test <- predict(xgb_CD8_TN_CN,CN_CD8_Tn_test)
xgb_CD8_TN_CN_predict_test <- data.frame(xgb_CD8_TN_CN_predict_test)
print(lm(xgb_CD8_TN_CN_predict_test$xgb_CD8_TN_CN_predict_test~CN_CD8_Tn_test$CN))
print(cor(xgb_CD8_TN_CN_predict_test$xgb_CD8_TN_CN_predict_test,CN_CD8_Tn_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
# gc()
# rm(xgb_CD8_TN_CN)
#
#
#
#
#
#
# ###__________________________________________________________________________________________###
# ###____________________________________________CD8_Tcm________________________________________###
# ###__________________________________________________________________________________________###
#
#
# ## prep data ##
# CN_CD8_Tcm <- read.delim("~/Analysis/data/protein/CD8_Tcm_CN_log10.csv",sep=";",dec = ",")
# CN_CD8_Tcm <- subset(CN_CD8_Tcm,CN_CD8_Tcm$CN>0)
# dim(CN_CD8_Tcm)
#
# CN_CD8_Tcm <- merge(CN_CD8_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_CD8_Tcm[is.na(CN_CD8_Tcm)]= 0
# print(dim(CN_CD8_Tcm)) # 7984 7126
#
#
# rownames(CN_CD8_Tcm) <- CN_CD8_Tcm$ID
# CN_CD8_Tcm$ID <- NULL
#
#
# ## Feature selection ##
# registerDoMC(6)
#
# CN_CD8_Tcm[is.na(CN_CD8_Tcm)]=0
# CN_CD8_Tcm <- CN_CD8_Tcm[, -nearZeroVar(CN_CD8_Tcm, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_CD8_Tcm)) # 7984 2759
#
# gc()
#
#
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_CD8_Tcm), 0.8 * nrow(CN_CD8_Tcm))
# CN_CD8_Tcm_train <- CN_CD8_Tcm[train_row,]
# CN_CD8_Tcm_test <- CN_CD8_Tcm[-train_row,]
#
# dim(CN_CD8_Tcm_train)
# dim(CN_CD8_Tcm_test)
#
# write.table(CN_CD8_Tcm_train,"~/Analysis/split_dataset/CN_CD8_Tcm_train",sep=";")
# write.table(CN_CD8_Tcm_test,"~/Analysis/split_dataset/CN_CD8_Tcm_test",sep=";")
#
CN_CD8_Tcm_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tcm_train",sep=";")
CN_CD8_Tcm_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tcm_test",sep=";")
#
#
#
# ## Training the xgb model ##
# print("modeling CD8 TN CN")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_CD8_Tcm_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                  max_depth = 6,
#                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                  ## The values below are default values in the sklearn-api.
#                                  eta = 0.05,
#                                  gamma=0,
#                                  min_child_weight = 0.9,
#                                  subsample = 1)
#
# start_time <- Sys.time()
# xgb_CD8_Tcm_CN <- train(CN~.,
#                        data=CN_CD8_Tcm_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_CD8_Tcm_CN,
#                        na.action = na.omit,
#                        nthread=72,
#                        verbose = TRUE)
#
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
#
# saveRDS(xgb_CD8_Tcm_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tcm_18_10_2021.RDS")
xgb_CD8_Tcm_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tcm_18_10_2021.RDS")
#
# print(xgb_CD8_Tcm_CN)
#
xgb_CD8_Tcm_CN_predict_test <- predict(xgb_CD8_Tcm_CN,CN_CD8_Tcm_test)
xgb_CD8_Tcm_CN_predict_test <- data.frame(xgb_CD8_Tcm_CN_predict_test)
print(lm(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test~CN_CD8_Tcm_test$CN))
print(cor(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test,CN_CD8_Tcm_test$CN, method = "pearson", use = "complete.obs")^2)

#plot(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test,CN_CD8_Tcm_test$CN)+abline(lm(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test~CN_CD8_Tcm_test$CN))
#
# gc()
# rm(xgb_CD8_Tcm_CN,CN_CD8_Tcm)
#
#
#
#
#
#
# ###__________________________________________________________________________________________###
# ###____________________________________________CD8_Tem________________________________________###
# ###__________________________________________________________________________________________###
#
#
# ## prep data ##
# CN_CD8_Tem <- read.delim("~/Analysis/data/protein/CD8_Tem_CN_log10.csv",sep=";",dec = ",")
# CN_CD8_Tem <- subset(CN_CD8_Tem,CN_CD8_Tem$CN>0)
# dim(CN_CD8_Tem)
#
# CN_CD8_Tem <- merge(CN_CD8_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_CD8_Tem[is.na(CN_CD8_Tem)]= 0
# print(dim(CN_CD8_Tem)) # 7984 7126
#
#
# rownames(CN_CD8_Tem) <- CN_CD8_Tem$ID
# CN_CD8_Tem$ID <- NULL
#
#
# ## Feature selection ##
# registerDoMC(6)
#
# CN_CD8_Tem[is.na(CN_CD8_Tem)]=0
# CN_CD8_Tem <- CN_CD8_Tem[, -nearZeroVar(CN_CD8_Tem, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_CD8_Tem)) # 7984 2759
#
# gc()
#
#
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_CD8_Tem), 0.8 * nrow(CN_CD8_Tem))
# CN_CD8_Tem_train <- CN_CD8_Tem[train_row,]
# CN_CD8_Tem_test <- CN_CD8_Tem[-train_row,]
#
# dim(CN_CD8_Tem_train)
# dim(CN_CD8_Tem_test)
#
# write.table(CN_CD8_Tem_train,"~/Analysis/split_dataset/CN_CD8_Tem_train",sep=";")
# write.table(CN_CD8_Tem_test,"~/Analysis/split_dataset/CN_CD8_Tem_test",sep=";")
#
CN_CD8_Tem_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tem_train",sep=";")
CN_CD8_Tem_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tem_test",sep=";")
#
#
#
# ## Training the xgb model ##
# print("modeling CD8 TN CN")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_CD8_Tem_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                   max_depth = 6,
#                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                   ## The values below are default values in the sklearn-api.
#                                   eta = 0.05,
#                                   gamma=0,
#                                   min_child_weight = 0.9,
#                                   subsample = 1)
#
# start_time <- Sys.time()
# xgb_CD8_Tem_CN <- train(CN~.,
#                         data=CN_CD8_Tem_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_CD8_Tem_CN,
#                         na.action = na.omit,
#                         nthread=72,
#                         verbose = TRUE)
#
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
#
# print(xgb_CD8_Tem_CN)
#
# saveRDS(xgb_CD8_Tem_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tem_18_10_2021.RDS")
xgb_CD8_Tem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tem_18_10_2021.RDS")
#
xgb_CD8_Tem_CN_predict_test <- predict(xgb_CD8_Tem_CN,CN_CD8_Tem_test)
xgb_CD8_Tem_CN_predict_test <- data.frame(xgb_CD8_Tem_CN_predict_test)
print(lm(xgb_CD8_Tem_CN_predict_test$xgb_CD8_Tem_CN_predict_test~CN_CD8_Tem_test$CN))
print(cor(xgb_CD8_Tem_CN_predict_test$xgb_CD8_Tem_CN_predict_test,CN_CD8_Tem_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
# gc()
# rm(xgb_CD8_Tem_CN,CN_CD8_Tem)
#
#
#
#
#
#
#
# ###__________________________________________________________________________________________###
# ###____________________________________________CD8_Teff________________________________________###
# ###__________________________________________________________________________________________###
#
#
# ## prep data ##
# CN_CD8_Teff <- read.delim("~/Analysis/data/protein/CD8_Teff_CN_log10.csv",sep=";",dec = ",")
# CN_CD8_Teff <- subset(CN_CD8_Teff,CN_CD8_Teff$CN>0)
# dim(CN_CD8_Teff)
#
# CN_CD8_Teff <- merge(CN_CD8_Teff,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_CD8_Teff[is.na(CN_CD8_Teff)]= 0
# print(dim(CN_CD8_Teff)) # 7984 7126
#
#
# rownames(CN_CD8_Teff) <- CN_CD8_Teff$ID
# CN_CD8_Teff$ID <- NULL
#
#
# ## Feature selection ##
# registerDoMC(6)
#
# CN_CD8_Teff[is.na(CN_CD8_Teff)]=0
# CN_CD8_Teff <- CN_CD8_Teff[, -nearZeroVar(CN_CD8_Teff, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_CD8_Teff)) # 7984 2759
#
# gc()
#
#
#
#
#
#
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_CD8_Teff), 0.8 * nrow(CN_CD8_Teff))
# CN_CD8_Teff_train <- CN_CD8_Teff[train_row,]
# CN_CD8_Teff_test <- CN_CD8_Teff[-train_row,]
#
# dim(CN_CD8_Teff_train)
# dim(CN_CD8_Teff_test)
#
# write.table(CN_CD8_Teff_train,"~/Analysis/split_dataset/CN_CD8_Teff_train",sep=";")
# write.table(CN_CD8_Teff_test,"~/Analysis/split_dataset/CN_CD8_Teff_test",sep=";")
#
CN_CD8_Teff_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Teff_train",sep=";")
CN_CD8_Teff_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Teff_test",sep=";")
#
#
#
# ## Training the xgb model ##
# print("modeling CD8 TN CN")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
#
# xgbGrid_CD8_Teff_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                   max_depth = 6,
#                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                   ## The values below are default values in the sklearn-api.
#                                   eta = 0.05,
#                                   gamma=0,
#                                   min_child_weight = 0.9,
#                                   subsample = 1)
#
# start_time <- Sys.time()
# xgb_CD8_Teff_CN <- train(CN~.,
#                         data=CN_CD8_Teff_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_CD8_Teff_CN,
#                         na.action = na.omit,
#                         nthread=72,
#                         verbose = TRUE)
#
#
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
#
# print(xgb_CD8_Teff_CN)
#
# saveRDS(xgb_CD8_Teff_CN, "~/Analysis/models/protein/xgb_CN_CD8_Teff_18_10_2021.RDS")
xgb_CD8_Teff_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Teff_18_10_2021.RDS")
xgb_CD8_Teff_CN_predict_test <- predict(xgb_CD8_Teff_CN,CN_CD8_Teff_test)
xgb_CD8_Teff_CN_predict_test <- data.frame(xgb_CD8_Teff_CN_predict_test)
print(lm(xgb_CD8_Teff_CN_predict_test$xgb_CD8_Teff_CN_predict_test~CN_CD8_Teff_test$CN))
print(cor(xgb_CD8_Teff_CN_predict_test$xgb_CD8_Teff_CN_predict_test,CN_CD8_Teff_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
# gc()
# rm(xgb_CD8_Teff_CN,CN_CD8_Teff)





###__________________________________________________________________________________________###
###____________________________________________CD4_Tn________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tn <- read.delim("~/Analysis/data/protein/CD4_Tn_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tn <- subset(CN_CD4_Tn,CN_CD4_Tn$CN>0)
dim(CN_CD4_Tn)

CN_CD4_Tn <- merge(CN_CD4_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tn[is.na(CN_CD4_Tn)]= 0
print(dim(CN_CD4_Tn)) # 7984 7126


rownames(CN_CD4_Tn) <- CN_CD4_Tn$ID
CN_CD4_Tn$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tn[is.na(CN_CD4_Tn)]=0
CN_CD4_Tn <- CN_CD4_Tn[, -nearZeroVar(CN_CD4_Tn, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD4_Tn)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD4_Tn), 0.8 * nrow(CN_CD4_Tn))
CN_CD4_Tn_train <- CN_CD4_Tn[train_row,]
CN_CD4_Tn_test <- CN_CD4_Tn[-train_row,]

dim(CN_CD4_Tn_train)
dim(CN_CD4_Tn_test)

write.table(CN_CD4_Tn_train,"~/Analysis/split_dataset/CN_CD4_Tn_train",sep=";")
write.table(CN_CD4_Tn_test,"~/Analysis/split_dataset/CN_CD4_Tn_test",sep=";")

# CN_CD4_Tn_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_train",sep=";")
# CN_CD4_Tn_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_test",sep=";")



## Training the xgb model ##
print("modeling CD4 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4_Tn_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_CD4_Tn_CN <- train(CN~.,
                        data=CN_CD4_Tn_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD4_Tn_CN,
                        na.action = na.omit,
                        nthread=94,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(xgb_CD4_Tn_CN)

saveRDS(xgb_CD4_Tn_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tn_18_10_2021.RDS")
xgb_CD4_Tn_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tn_18_10_2021.RDS")

xgb_CD4_Tn_CN_predict_test <- predict(xgb_CD4_Tn_CN,CN_CD4_Tn_test)
xgb_CD4_Tn_CN_predict_test <- data.frame(xgb_CD4_Tn_CN_predict_test)
print(lm(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test~CN_CD4_Tn_test$CN))

print(cor(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN)


gc()
rm(xgb_CD4_Tn_CN,CN_CD4_Tn)









###__________________________________________________________________________________________###
###____________________________________________CD4_Tcm________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tcm <- read.delim("~/Analysis/data/protein/CD4_Tcm_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tcm <- subset(CN_CD4_Tcm,CN_CD4_Tcm$CN>0)
dim(CN_CD4_Tcm)

CN_CD4_Tcm <- merge(CN_CD4_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tcm[is.na(CN_CD4_Tcm)]= 0
print(dim(CN_CD4_Tcm)) # 7984 7126


rownames(CN_CD4_Tcm) <- CN_CD4_Tcm$ID
CN_CD4_Tcm$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tcm[is.na(CN_CD4_Tcm)]=0
CN_CD4_Tcm <- CN_CD4_Tcm[, -nearZeroVar(CN_CD4_Tcm, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD4_Tcm)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD4_Tcm), 0.8 * nrow(CN_CD4_Tcm))
CN_CD4_Tcm_train <- CN_CD4_Tcm[train_row,]
CN_CD4_Tcm_test <- CN_CD4_Tcm[-train_row,]

dim(CN_CD4_Tcm_train)
dim(CN_CD4_Tcm_test)

write.table(CN_CD4_Tcm_train,"~/Analysis/split_dataset/CN_CD4_Tcm_train",sep=";")
write.table(CN_CD4_Tcm_test,"~/Analysis/split_dataset/CN_CD4_Tcm_test",sep=";")

# CN_CD4_Tcm_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tcm_train",sep=";")
# CN_CD4_Tcm_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tcm_test",sep=";")



## Training the xgb model ##
print("modeling CD4 Tcm CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4_Tcm_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                 max_depth = 6,
                                 colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                 ## The values below are default values in the sklearn-api.
                                 eta = 0.05,
                                 gamma=0,
                                 min_child_weight = 0.9,
                                 subsample = 1)

start_time <- Sys.time()
xgb_CD4_Tcm_CN <- train(CN~.,
                       data=CN_CD4_Tcm_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_CD4_Tcm_CN,
                       na.action = na.omit,
                       nthread=72,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(xgb_CD4_Tcm_CN)

saveRDS(xgb_CD4_Tcm_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tcm_18_10_2021.RDS")
xgb_CD4_Tcm_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tcm_18_10_2021.RDS")

xgb_CD4_Tcm_CN_predict_test <- predict(xgb_CD4_Tcm_CN,CN_CD4_Tcm_test)
xgb_CD4_Tcm_CN_predict_test <- data.frame(xgb_CD4_Tcm_CN_predict_test)
print(lm(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test~CN_CD4_Tcm_test$CN))

print(cor(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test,CN_CD4_Tcm_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test,CN_CD4_Tcm_test$CN)


gc()
rm(xgb_CD4_Tcm_CN,CN_CD4_Tcm)









###__________________________________________________________________________________________###
###____________________________________________CD4_Tem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tem <- read.delim("~/Analysis/data/protein/CD4_Tem_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tem <- subset(CN_CD4_Tem,CN_CD4_Tem$CN>0)
dim(CN_CD4_Tem)

CN_CD4_Tem <- merge(CN_CD4_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tem[is.na(CN_CD4_Tem)]= 0
print(dim(CN_CD4_Tem)) # 7984 7126


rownames(CN_CD4_Tem) <- CN_CD4_Tem$ID
CN_CD4_Tem$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tem[is.na(CN_CD4_Tem)]=0
CN_CD4_Tem <- CN_CD4_Tem[, -nearZeroVar(CN_CD4_Tem, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD4_Tem)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD4_Tem), 0.8 * nrow(CN_CD4_Tem))
CN_CD4_Tem_train <- CN_CD4_Tem[train_row,]
CN_CD4_Tem_test <- CN_CD4_Tem[-train_row,]

dim(CN_CD4_Tem_train)
dim(CN_CD4_Tem_test)

write.table(CN_CD4_Tem_train,"~/Analysis/split_dataset/CN_CD4_Tem_train",sep=";")
write.table(CN_CD4_Tem_test,"~/Analysis/split_dataset/CN_CD4_Tem_test",sep=";")

CN_CD4_Tem_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_train",sep=";")
CN_CD4_Tem_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_test",sep=";")



## Training the xgb model ##
print("modeling CD4 Tem CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4_Tem_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                 max_depth = 6,
                                 colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                 ## The values below are default values in the sklearn-api.
                                 eta = 0.05,
                                 gamma=0,
                                 min_child_weight = 0.9,
                                 subsample = 1)

start_time <- Sys.time()
xgb_CD4_Tem_CN <- train(CN~.,
                       data=CN_CD4_Tem_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_CD4_Tem_CN,
                       na.action = na.omit,
                       nthread=72,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(xgb_CD4_Tem_CN)

saveRDS(xgb_CD4_Tem_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tem_18_10_2021.RDS")
xgb_CD4_Tem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tem_18_10_2021.RDS")

xgb_CD4_Tem_CN_predict_test <- predict(xgb_CD4_Tem_CN,CN_CD4_Tem_test)
xgb_CD4_Tem_CN_predict_test <- data.frame(xgb_CD4_Tem_CN_predict_test)
print(lm(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test~CN_CD4_Tem_test$CN))

print(cor(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test,CN_CD4_Tem_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test,CN_CD4_Tem_test$CN)


gc()
rm(xgb_CD4_Tem_CN,CN_CD4_Tem)





