### protein models including mRNA measurements ###

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
###____________________________________________BiomaRt________________________________________###
###__________________________________________________________________________________________###



## Biomart ##
# ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")

ensembl <- useMart("ensembl","hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")

# ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)
tx2gene <- subset(tx2gene,tx2gene$transcript_biotype == "protein_coding")
tx2gene$ccds <- NULL
tx2gene$ensembl_transcript_id_version <- NULL

attr_ens <- listAttributes(ensembl)
gene2PtID <- getBM(attributes=c("ensembl_gene_id","uniprotswissprot"), mart = ensembl)
#gene2PtID$uniprot_isoform <- mapply(strsplit(as.character(gene2PtID$uniprot_isoform),"-"),FUN=function(x){(as.character(x)[1])})

gene2PtID <- gene2PtID[gene2PtID$uniprotswissprot>1,]


###__________________________________________________________________________________________###
###____________________________________________ex vivo________________________________________###
###__________________________________________________________________________________________###


## parameters ##
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")



###__________________________________________________________________________________________###
###____________________________________________HEK293________________________________________###
###__________________________________________________________________________________________###

# 
# ## prep data ##
# CN_HEK293 <- read.delim("~/Analysis/data/HEK293_CN_log10.csv",sep=";",dec = ",")
# # CN_HEK293 <- subset(CN_HEK293,CN_HEK293$CN>0)
# dim(CN_HEK293)
# 
# 
# CN_HEK293 <- merge(CN_HEK293,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_HEK293[is.na(CN_HEK293)]= 0
# print(dim(CN_HEK293)) # 7984 7126
# 
# CN_HEK293_dups <- CN_HEK293[duplicated(CN_HEK293$ID),]
# 
# rownames(CN_HEK293) <- CN_HEK293$ID
# CN_HEK293$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_HEK293[is.na(CN_HEK293)]=0
# CN_HEK293 <- CN_HEK293[, -nearZeroVar(CN_HEK293, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_HEK293)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_HEK293), 0.8 * nrow(CN_HEK293))
# CN_HEK293_train <- CN_HEK293[train_row,]
# CN_HEK293_test <- CN_HEK293[-train_row,]
# 
# dim(CN_HEK293_train)
# dim(CN_HEK293_test)
# 
# write.table(CN_HEK293_train,"~/Analysis/split_dataset/CN_HEK293_train",sep=";")
# write.table(CN_HEK293_test,"~/Analysis/split_dataset/CN_HEK293_test",sep=";")
# 
# # CN_HEK293_train <- read.delim("~/Analysis/split_dataset/CN_HEK293_train",sep=";")
# # CN_HEK293_test <- read.delim("~/Analysis/split_dataset/CN_HEK293_test",sep=";")
# #
# #
# #
# # ## Training the xgb model ##
# print("modeling HEK293")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_HEK293_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                               max_depth = 6,
#                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                               ## The values below are default values in the sklearn-api.
#                               eta = 0.05,
#                               gamma=0,
#                               min_child_weight = 0.9,
#                               subsample = 1)
# 
# start_time <- Sys.time()
# xgb_HEK293_CN <- train(CN~.,
#                         data=CN_HEK293_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_HEK293_CN,
#                         na.action = na.omit,
#                         nthread=80,
#                         verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_HEK293_CN, "~/Analysis/models/protein/xgb_CN_HEK293_07_06_2022.RDS")
# # xgb_HEK293_CN <- readRDS("~/Analysis/models/protein/xgb_CN_HEK293_07_06_2022.RDS")
# #
# # print(xgb_HEK293_CN)
# 
# xgb_HEK293_CN_predict_test <- predict(xgb_HEK293_CN,CN_HEK293_test)
# xgb_HEK293_CN_predict_test <- data.frame(xgb_HEK293_CN_predict_test)
# # print(lm(xgb_HEK293_CN_predict_test$xgb_HEK293_CN_predict_test~CN_HEK293_test$CN))
# print(cor(xgb_HEK293_CN_predict_test$xgb_HEK293_CN_predict_test,CN_HEK293_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(CN_HEK293_train$TPM,CN_HEK293_train$CN)
# # print(cor(CN_HEK293_train$TPM,CN_HEK293_train$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(xgb_HEK293_CN_predict_test$xgb_HEK293_CN_predict_test,CN_HEK293_test$CN)
# 
# gc()
# rm(xgb_HEK293_CN)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###__________________________________________________________________________________________###
# ###____________________________________________HeLa________________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# CN_HeLa <- read.delim("~/Analysis/data/HeLa_CN_log10.csv",sep=";",dec = ",")
# # CN_HeLa <- subset(CN_HeLa,CN_HeLa$CN>0)
# dim(CN_HeLa)
# 
# 
# CN_HeLa <- merge(CN_HeLa,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_HeLa[is.na(CN_HeLa)]= 0
# print(dim(CN_HeLa)) # 7984 7126
# 
# CN_HeLa_dups <- CN_HeLa[duplicated(CN_HeLa$ID),]
# 
# rownames(CN_HeLa) <- CN_HeLa$ID
# CN_HeLa$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_HeLa[is.na(CN_HeLa)]=0
# CN_HeLa <- CN_HeLa[, -nearZeroVar(CN_HeLa, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_HeLa)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_HeLa), 0.8 * nrow(CN_HeLa))
# CN_HeLa_train <- CN_HeLa[train_row,]
# CN_HeLa_test <- CN_HeLa[-train_row,]
# 
# dim(CN_HeLa_train)
# dim(CN_HeLa_test)
# 
# write.table(CN_HeLa_train,"~/Analysis/split_dataset/CN_HeLa_train",sep=";")
# write.table(CN_HeLa_test,"~/Analysis/split_dataset/CN_HeLa_test",sep=";")
# 
# # CN_HeLa_train <- read.delim("~/Analysis/split_dataset/CN_HeLa_train",sep=";")
# # CN_HeLa_test <- read.delim("~/Analysis/split_dataset/CN_HeLa_test",sep=";")
# #
# #
# #
# # ## Training the xgb model ##
# print("modeling HeLa")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_HeLa_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                  max_depth = 6,
#                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                  ## The values below are default values in the sklearn-api.
#                                  eta = 0.05,
#                                  gamma=0,
#                                  min_child_weight = 0.9,
#                                  subsample = 1)
# 
# start_time <- Sys.time()
# xgb_HeLa_CN <- train(CN~.,
#                        data=CN_HeLa_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_HeLa_CN,
#                        na.action = na.omit,
#                        nthread=80,
#                        verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_HeLa_CN, "~/Analysis/models/protein/xgb_CN_HeLa_07_06_2022.RDS")
# # xgb_HeLa_CN <- readRDS("~/Analysis/models/protein/xgb_CN_HeLa_07_06_2022.RDS")
# #
# # print(xgb_HeLa_CN)
# 
# xgb_HeLa_CN_predict_test <- predict(xgb_HeLa_CN,CN_HeLa_test)
# xgb_HeLa_CN_predict_test <- data.frame(xgb_HeLa_CN_predict_test)
# # print(lm(xgb_HeLa_CN_predict_test$xgb_HeLa_CN_predict_test~CN_HeLa_test$CN))
# print(cor(xgb_HeLa_CN_predict_test$xgb_HeLa_CN_predict_test,CN_HeLa_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(CN_HeLa_train$TPM,CN_HeLa_train$CN)
# # print(cor(CN_HeLa_train$TPM,CN_HeLa_train$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(xgb_HeLa_CN_predict_test$xgb_HeLa_CN_predict_test,CN_HeLa_test$CN)
# 
# gc()
# rm(xgb_HeLa_CN)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###__________________________________________________________________________________________###
# ###____________________________________________K562________________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# CN_K562 <- read.delim("~/Analysis/data/K562_CN_log10.csv",sep=";",dec = ",")
# # CN_K562 <- subset(CN_K562,CN_K562$CN>0)
# dim(CN_K562)
# 
# 
# CN_K562 <- merge(CN_K562,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_K562[is.na(CN_K562)]= 0
# print(dim(CN_K562)) # 7984 7126
# 
# CN_K562_dups <- CN_K562[duplicated(CN_K562$ID),]
# 
# rownames(CN_K562) <- CN_K562$ID
# CN_K562$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_K562[is.na(CN_K562)]=0
# CN_K562 <- CN_K562[, -nearZeroVar(CN_K562, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_K562)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_K562), 0.8 * nrow(CN_K562))
# CN_K562_train <- CN_K562[train_row,]
# CN_K562_test <- CN_K562[-train_row,]
# 
# dim(CN_K562_train)
# dim(CN_K562_test)
# 
# write.table(CN_K562_train,"~/Analysis/split_dataset/CN_K562_train",sep=";")
# write.table(CN_K562_test,"~/Analysis/split_dataset/CN_K562_test",sep=";")
# 
# # CN_K562_train <- read.delim("~/Analysis/split_dataset/CN_K562_train",sep=";")
# # CN_K562_test <- read.delim("~/Analysis/split_dataset/CN_K562_test",sep=";")
# #
# #
# #
# # ## Training the xgb model ##
# print("modeling K562")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_K562_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                  max_depth = 6,
#                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                  ## The values below are default values in the sklearn-api.
#                                  eta = 0.05,
#                                  gamma=0,
#                                  min_child_weight = 0.9,
#                                  subsample = 1)
# 
# start_time <- Sys.time()
# xgb_K562_CN <- train(CN~.,
#                        data=CN_K562_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_K562_CN,
#                        na.action = na.omit,
#                        nthread=80,
#                        verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_K562_CN, "~/Analysis/models/protein/xgb_CN_K562_07_06_2022.RDS")
# # xgb_K562_CN <- readRDS("~/Analysis/models/protein/xgb_CN_K562_07_06_2022.RDS")
# #
# # print(xgb_K562_CN)
# 
# xgb_K562_CN_predict_test <- predict(xgb_K562_CN,CN_K562_test)
# xgb_K562_CN_predict_test <- data.frame(xgb_K562_CN_predict_test)
# # print(lm(xgb_K562_CN_predict_test$xgb_K562_CN_predict_test~CN_K562_test$CN))
# print(cor(xgb_K562_CN_predict_test$xgb_K562_CN_predict_test,CN_K562_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(CN_K562_train$TPM,CN_K562_train$CN)
# # print(cor(CN_K562_train$TPM,CN_K562_train$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(xgb_K562_CN_predict_test$xgb_K562_CN_predict_test,CN_K562_test$CN)
# 
# gc()
# rm(xgb_K562_CN)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##__________________________________________________________________________________________________________##
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###__________________________________________________________________________________________###
# ###____________________________________________HEK293________________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# CN_HEK293 <- read.delim("~/Analysis/data/HEK293_CN_log10.csv",sep=";",dec = ",")
# CN_HEK293 <- subset(CN_HEK293,CN_HEK293$CN>0)
# dim(CN_HEK293)
# 
# 
# TPM_HEK293 <- read.delim("~/Analysis/data/HEK293T_TPM_log10.csv",sep=";",dec=",")
# TPM_HEK293$TPM <- 10^TPM_HEK293$TPM
# TPM_HEK293 <- merge(TPM_HEK293,tx2gene, by.x="ID", by.y="ensembl_transcript_id")
# 
# 
# registerDoMC(10)
# 
# TPM_HEK293 <- merge(TPM_HEK293,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
# TPM_HEK293 <- ddply(TPM_HEK293,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
# TPM_HEK293$TPM <- log10(TPM_HEK293$TPM)
# #dim(TPM_HEK293[duplicated(TPM_HEK293$uniprotswissprot),])
# 
# 
# 
# CN_HEK293 <- merge(CN_HEK293,TPM_HEK293,by.x="ID",by.y="uniprotswissprot")
# CN_HEK293$ensembl_gene_id <- NULL
# CN_HEK293 <- merge(CN_HEK293,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_HEK293[is.na(CN_HEK293)]= 0
# print(dim(CN_HEK293)) # 7984 7126
# 
# CN_HEK293_dups <- CN_HEK293[duplicated(CN_HEK293$ID),]
# 
# rownames(CN_HEK293) <- CN_HEK293$ID
# CN_HEK293$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_HEK293[is.na(CN_HEK293)]=0
# CN_HEK293 <- CN_HEK293[, -nearZeroVar(CN_HEK293, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_HEK293)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_HEK293), 0.8 * nrow(CN_HEK293))
# CN_HEK293_train <- CN_HEK293[train_row,]
# CN_HEK293_test <- CN_HEK293[-train_row,]
# 
# dim(CN_HEK293_train)
# dim(CN_HEK293_test)
# 
# write.table(CN_HEK293_train,"~/Analysis/split_dataset/CN_HEK293_with_RNA_train",sep=";")
# write.table(CN_HEK293_test,"~/Analysis/split_dataset/CN_HEK293_with_RNA_test",sep=";")
# 
# # CN_HEK293_train <- read.delim("~/Analysis/split_dataset/CN_HEK293_with_RNA_train",sep=";")
# # CN_HEK293_test <- read.delim("~/Analysis/split_dataset/CN_HEK293_with_RNA_test",sep=";")
# #
# #
# #
# # ## Training the xgb model ##
# print("modeling HEK293 with RNA")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_HEK293_CN_w_RNA <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                               max_depth = 6,
#                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                               ## The values below are default values in the sklearn-api.
#                               eta = 0.05,
#                               gamma=0,
#                               min_child_weight = 0.9,
#                               subsample = 1)
# 
# start_time <- Sys.time()
# xgb_HEK293_CN_w_RNA <- train(CN~.,
#                         data=CN_HEK293_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_HEK293_CN_w_RNA,
#                         na.action = na.omit,
#                         nthread=80,
#                         verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_HEK293_CN_w_RNA, "~/Analysis/models/protein/xgb_CN_HEK293_with_RNA_07_06_2022.RDS")
# # xgb_HEK293_CN_w_RNA <- readRDS("~/Analysis/models/protein/xgb_CN_HEK293_with_RNA_07_06_2022.RDS")
# #
# # print(xgb_HEK293_CN_w_RNA)
# 
# xgb_HEK293_CN_w_RNA_predict_test <- predict(xgb_HEK293_CN_w_RNA,CN_HEK293_test)
# xgb_HEK293_CN_w_RNA_predict_test <- data.frame(xgb_HEK293_CN_w_RNA_predict_test)
# # print(lm(xgb_HEK293_CN_w_RNA_predict_test$xgb_HEK293_CN_w_RNA_predict_test~CN_HEK293_test$CN))
# print(cor(xgb_HEK293_CN_w_RNA_predict_test$xgb_HEK293_CN_w_RNA_predict_test,CN_HEK293_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(CN_HEK293_with_RNA_train$TPM,CN_HEK293_with_RNA_train$CN)
# print(cor(CN_HEK293_train$TPM,CN_HEK293_train$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(xgb_HEK293_CN_w_RNA_predict_test$xgb_HEK293_CN_w_RNA_predict_test,CN_HEK293_test$CN)
# 
# gc()
# rm(xgb_HEK293_CN_w_RNA)
# 
# 
# 
# 
# 
# 
# 
# 
# ###__________________________________________________________________________________________###
# ###____________________________________________HeLa________________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# CN_HeLa <- read.delim("~/Analysis/data/HeLa_CN_log10.csv",sep=";",dec = ",")
# CN_HeLa <- subset(CN_HeLa,CN_HeLa$CN>0)
# dim(CN_HeLa)
# 
# 
# TPM_HeLa <- read.delim("~/Analysis/data/HeLa_TPM_log10.csv",sep=";",dec=",")
# TPM_HeLa$TPM <- 10^TPM_HeLa$TPM
# TPM_HeLa <- merge(TPM_HeLa,tx2gene, by.x="ID", by.y="ensembl_transcript_id")
# 
# 
# registerDoMC(10)
# 
# TPM_HeLa <- merge(TPM_HeLa,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
# TPM_HeLa <- ddply(TPM_HeLa,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
# TPM_HeLa$TPM <- log10(TPM_HeLa$TPM)
# #dim(TPM_HeLa[duplicated(TPM_HeLa$uniprotswissprot),])
# 
# 
# 
# CN_HeLa <- merge(CN_HeLa,TPM_HeLa,by.x="ID",by.y="uniprotswissprot")
# CN_HeLa$ensembl_gene_id <- NULL
# CN_HeLa <- merge(CN_HeLa,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_HeLa[is.na(CN_HeLa)]= 0
# print(dim(CN_HeLa)) # 7984 7126
# 
# CN_HeLa_dups <- CN_HeLa[duplicated(CN_HeLa$ID),]
# 
# rownames(CN_HeLa) <- CN_HeLa$ID
# CN_HeLa$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_HeLa[is.na(CN_HeLa)]=0
# CN_HeLa <- CN_HeLa[, -nearZeroVar(CN_HeLa, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_HeLa)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_HeLa), 0.8 * nrow(CN_HeLa))
# CN_HeLa_train <- CN_HeLa[train_row,]
# CN_HeLa_test <- CN_HeLa[-train_row,]
# 
# dim(CN_HeLa_train)
# dim(CN_HeLa_test)
# 
# write.table(CN_HeLa_train,"~/Analysis/split_dataset/CN_HeLa_with_RNA_train",sep=";")
# write.table(CN_HeLa_test,"~/Analysis/split_dataset/CN_HeLa_with_RNA_test",sep=";")
# 
# # CN_HeLa_train <- read.delim("~/Analysis/split_dataset/CN_HeLa_with_RNA_train",sep=";")
# # CN_HeLa_test <- read.delim("~/Analysis/split_dataset/CN_HeLa_with_RNA_test",sep=";")
# #
# #
# #
# # ## Training the xgb model ##
# print("modeling HeLa with RNA")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_HeLa_CN_w_RNA <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                                        max_depth = 6,
#                                        colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                                        ## The values below are default values in the sklearn-api.
#                                        eta = 0.05,
#                                        gamma=0,
#                                        min_child_weight = 0.9,
#                                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_HeLa_CN_w_RNA <- train(CN~.,
#                        data=CN_HeLa_train,
#                        method="xgbTree",
#                        trControl=control,
#                        #metric="Rsquared",
#                        tuneGrid= xgbGrid_HeLa_CN_w_RNA,
#                        na.action = na.omit,
#                        nthread=80,
#                        verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_HeLa_CN_w_RNA, "~/Analysis/models/protein/xgb_CN_HeLa_with_RNA_07_06_2022.RDS")
# # xgb_HeLa_CN_w_RNA <- readRDS("~/Analysis/models/protein/xgb_CN_HeLa_with_RNA_07_06_2022.RDS")
# #
# # print(xgb_HeLa_CN_w_RNA)
# 
xgb_HeLa_CN_w_RNA_predict_test <- predict(xgb_HeLa_CN_w_RNA,CN_HeLa_test)
xgb_HeLa_CN_w_RNA_predict_test <- data.frame(xgb_HeLa_CN_w_RNA_predict_test)
# print(lm(xgb_HeLa_CN_w_RNA_predict_test$xgb_HeLa_CN_w_RNA_predict_test~CN_HeLa_test$CN))
print(cor(xgb_HeLa_CN_w_RNA_predict_test$xgb_HeLa_CN_w_RNA_predict_test,CN_HeLa_test$CN, method = "pearson", use = "complete.obs")^2)

plot(CN_HeLa_train$TPM,CN_HeLa_train$CN)
print(cor(CN_HeLa_train$TPM,CN_HeLa_train$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_HeLa_CN_w_RNA_predict_test$xgb_HeLa_CN_w_RNA_predict_test,CN_HeLa_test$CN,xlim = c(2,8),ylim = c(2,8),asp = 1)

# gc()
# rm(xgb_HeLa_CN_w_RNA)
# 
# 
# 
# 
# 
# 
# 
# 
# 





###__________________________________________________________________________________________###
###____________________________________________K562________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_K562 <- read.delim("~/Analysis/data/K562_CN_log10.csv",sep=";",dec = ",")
CN_K562 <- subset(CN_K562,CN_K562$CN>0)
dim(CN_K562)


TPM_K562 <- read.delim("~/Analysis/data/K562_TPM_log10.csv",sep=";",dec=",")
TPM_K562$TPM <- 10^TPM_K562$TPM
TPM_K562 <- merge(TPM_K562,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_K562 <- merge(TPM_K562,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_K562 <- ddply(TPM_K562,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_K562$TPM <- log10(TPM_K562$TPM)
#dim(TPM_K562[duplicated(TPM_K562$uniprotswissprot),])


CN_K562 <- merge(CN_K562,TPM_K562,by.x="ID",by.y="uniprotswissprot")
CN_K562$ensembl_gene_id <- NULL
CN_K562 <- merge(CN_K562,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_K562[is.na(CN_K562)]= 0
print(dim(CN_K562)) # 7984 7126

CN_K562_dups <- CN_K562[duplicated(CN_K562$ID),]

rownames(CN_K562) <- CN_K562$ID
CN_K562$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_K562[is.na(CN_K562)]=0
CN_K562 <- CN_K562[, -nearZeroVar(CN_K562, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_K562)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_K562), 0.8 * nrow(CN_K562))
CN_K562_train <- CN_K562[train_row,]
CN_K562_test <- CN_K562[-train_row,]

dim(CN_K562_train)
dim(CN_K562_test)

write.table(CN_K562_train,"~/Analysis/split_dataset/CN_K562_with_RNA_train",sep=";")
write.table(CN_K562_test,"~/Analysis/split_dataset/CN_K562_with_RNA_test",sep=";")

# CN_K562_train <- read.delim("~/Analysis/split_dataset/CN_K562_with_RNA_train",sep=";")
# CN_K562_test <- read.delim("~/Analysis/split_dataset/CN_K562_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling K562 with RNA")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_K562_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                       max_depth = 6,
                                       colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                       ## The values below are default values in the sklearn-api.
                                       eta = 0.05,
                                       gamma=0,
                                       min_child_weight = 0.9,
                                       subsample = 1)

start_time <- Sys.time()
xgb_K562_CN_w_RNA <- train(CN~.,
                       data=CN_K562_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_K562_CN,
                       na.action = na.omit,
                       nthread=90,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_K562_CN_w_RNA, "~/Analysis/models/protein/xgb_CN_K562_with_RNA_07_06_2022.RDS")
# xgb_K562_CN_w_RNA <- readRDS("~/Analysis/models/protein/xgb_CN_K562_with_RNA_07_06_2022.RDS")
#
# print(xgb_K562_CN_w_RNA)

xgb_K562_CN_w_RNA_predict_test <- predict(xgb_K562_CN_w_RNA,CN_K562_test)
xgb_K562_CN_w_RNA_predict_test <- data.frame(xgb_K562_CN_w_RNA_predict_test)
# print(lm(xgb_K562_CN_w_RNA_predict_test$xgb_K562_CN_w_RNA_predict_test~CN_K562_test$CN))
print(cor(xgb_K562_CN_w_RNA_predict_test$xgb_K562_CN_w_RNA_predict_test,CN_K562_test$CN, method = "pearson", use = "complete.obs")^2)

print(cor(CN_K562_train$TPM,CN_K562_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_K562_CN_w_RNA_predict_test$xgb_K562_CN_w_RNA_predict_test,CN_K562_test$CN)

gc()
rm(xgb_K562_CN_w_RNA)





