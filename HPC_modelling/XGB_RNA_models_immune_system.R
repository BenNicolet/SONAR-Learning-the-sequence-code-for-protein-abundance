library(plyr)
library(dplyr)
library(doMC)
library(biomaRt)
library(tidyverse)
library(caret)
library(e1071)
library(xgboost)




###__________________________________________________________________________________________###
###_______________________________________ immune subsets ___________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")

## Biomart ##
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)
#
#
# ## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###_________________________________________ B naive ________________________________________###
###__________________________________________________________________________________________###
# 
# ## prep data ##
# print("preping data for Tn 0h")
# TPM_Bnaive <- read.delim("~/Analysis/data/immune_cells_RNA/B_naive_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Bnaive <- merge(TPM_Bnaive,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Bnaive_param <- ddply(TPM_Bnaive,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Bnaive_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Bnaive$TPM <- 10^(TPM_Bnaive$TPM)
# TPM_Bnaive <- data.frame("ID"=TPM_Bnaive$ensembl_gene_id ,"TPM"=TPM_Bnaive$TPM)
# TPM_Bnaive <- ddply(TPM_Bnaive,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Bnaive <- merge(TPM_Bnaive,TPM_Bnaive_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Bnaive$TPM <- log10(TPM_Bnaive$TPM)
# 
# TPM_Bnaive[is.na(TPM_Bnaive)]= 0
# TPM_Bnaive <- subset(TPM_Bnaive, TPM_Bnaive$TPM>=-1)
# print(dim(TPM_Bnaive)) # 13729  2682
# 
# rownames(TPM_Bnaive) <- TPM_Bnaive$ID
# TPM_Bnaive$ID <- NULL
# TPM_Bnaive$ensembl_gene_id <- NULL
# TPM_Bnaive$gene_biotype <- NULL
# TPM_Bnaive$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Bnaive[is.na(TPM_Bnaive)]=0
# TPM_Bnaive <- TPM_Bnaive[, -caret::nearZeroVar(TPM_Bnaive, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Bnaive)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Bnaive), 0.8 * nrow(TPM_Bnaive))
# TPM_Bnaive_train <- TPM_Bnaive[train_row,]
# TPM_Bnaive_test <- TPM_Bnaive[-train_row,]
# 
# dim(TPM_Bnaive_train)
# dim(TPM_Bnaive_test)
# 
# write.table(TPM_Bnaive_train,"~/Analysis/split_dataset/TPM_Bnaive_train",sep=";")
# write.table(TPM_Bnaive_test,"~/Analysis/split_dataset/TPM_Bnaive_test",sep=";")
# 
# # TPM_Bnaive_train <- read.delim("~/Analysis/split_dataset/TPM_Bnaive_train",sep=";")
# # TPM_Bnaive_test <- read.delim("~/Analysis/split_dataset/TPM_Bnaive_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                              max_depth = 6,
#                              colsample_bytree = 0.4,
#                              eta = 0.05,
#                              gamma=0,
#                              min_child_weight = 0.9,
#                              subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Bnaive <- train(TPM~.,
#                         data=TPM_Bnaive_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Bnaive)
# 
# saveRDS(xgb_TPM_Bnaive, "~/Analysis/models/xgb_TPM_Bnaive_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_Bnaive_03_05_2022.RDS")
# 
# xgb_TPM_Bnaive_predict_test <- predict(xgb_TPM_Bnaive,TPM_Bnaive_test)
# xgb_TPM_Bnaive_predict_test <- data.frame(xgb_TPM_Bnaive_predict_test)
# print(lm(xgb_TPM_Bnaive_predict_test$xgb_TPM_Bnaive_predict_test~TPM_Bnaive_test$TPM))
# print(cor(xgb_TPM_Bnaive_predict_test$xgb_TPM_Bnaive_predict_test,TPM_Bnaive_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_Bnaive,TPM_Bnaive_param)
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
# ###__________________________________________________________________________________________###
# ###_________________________________________ B NS mem _______________________________________###
# ###__________________________________________________________________________________________###
# 
# ## prep data ##
# print("preping data for Tn 0h")
# TPM_B_NS_mem <- read.delim("~/Analysis/data/immune_cells_RNA/B_NS_mem_TPM_log10.csv",sep = ";",dec = ",")
# TPM_B_NS_mem <- merge(TPM_B_NS_mem,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_B_NS_mem_param <- ddply(TPM_B_NS_mem,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_B_NS_mem_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_B_NS_mem$TPM <- 10^(TPM_B_NS_mem$TPM)
# TPM_B_NS_mem <- data.frame("ID"=TPM_B_NS_mem$ensembl_gene_id ,"TPM"=TPM_B_NS_mem$TPM)
# TPM_B_NS_mem <- ddply(TPM_B_NS_mem,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_B_NS_mem <- merge(TPM_B_NS_mem,TPM_B_NS_mem_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_B_NS_mem$TPM <- log10(TPM_B_NS_mem$TPM)
# 
# TPM_B_NS_mem[is.na(TPM_B_NS_mem)]= 0
# TPM_B_NS_mem <- subset(TPM_B_NS_mem, TPM_B_NS_mem$TPM>=-1)
# print(dim(TPM_B_NS_mem)) # 13729  2682
# 
# rownames(TPM_B_NS_mem) <- TPM_B_NS_mem$ID
# TPM_B_NS_mem$ID <- NULL
# TPM_B_NS_mem$ensembl_gene_id <- NULL
# TPM_B_NS_mem$gene_biotype <- NULL
# TPM_B_NS_mem$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_B_NS_mem[is.na(TPM_B_NS_mem)]=0
# TPM_B_NS_mem <- TPM_B_NS_mem[, -caret::nearZeroVar(TPM_B_NS_mem, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_B_NS_mem)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_B_NS_mem), 0.8 * nrow(TPM_B_NS_mem))
# TPM_B_NS_mem_train <- TPM_B_NS_mem[train_row,]
# TPM_B_NS_mem_test <- TPM_B_NS_mem[-train_row,]
# 
# dim(TPM_B_NS_mem_train)
# dim(TPM_B_NS_mem_test)
# 
# write.table(TPM_B_NS_mem_train,"~/Analysis/split_dataset/TPM_B_NS_mem_train",sep=";")
# write.table(TPM_B_NS_mem_test,"~/Analysis/split_dataset/TPM_B_NS_mem_test",sep=";")
# 
# # TPM_B_NS_mem_train <- read.delim("~/Analysis/split_dataset/TPM_B_NS_mem_train",sep=";")
# # TPM_B_NS_mem_test <- read.delim("~/Analysis/split_dataset/TPM_B_NS_mem_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_B_NS_mem <- train(TPM~.,
#                         data=TPM_B_NS_mem_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_B_NS_mem)
# 
# saveRDS(xgb_TPM_B_NS_mem, "~/Analysis/models/xgb_TPM_B_NS_mem_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_B_NS_mem_03_05_2022.RDS")
# 
# xgb_TPM_B_NS_mem_predict_test <- predict(xgb_TPM_B_NS_mem,TPM_B_NS_mem_test)
# xgb_TPM_B_NS_mem_predict_test <- data.frame(xgb_TPM_B_NS_mem_predict_test)
# print(lm(xgb_TPM_B_NS_mem_predict_test$xgb_TPM_B_NS_mem_predict_test~TPM_B_NS_mem_test$TPM))
# print(cor(xgb_TPM_B_NS_mem_predict_test$xgb_TPM_B_NS_mem_predict_test,TPM_B_NS_mem_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_B_NS_mem,TPM_B_NS_mem_param)







###__________________________________________________________________________________________###
###_________________________________________ B plasma ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Tn 0h")
# TPM_B_plasma <- read.delim("~/Analysis/data/immune_cells_RNA/B_plasma_TPM_log10.csv",sep = ";",dec = ",")
# TPM_B_plasma <- merge(TPM_B_plasma,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_B_plasma_param <- ddply(TPM_B_plasma,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_B_plasma_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_B_plasma$TPM <- 10^(TPM_B_plasma$TPM)
# TPM_B_plasma <- data.frame("ID"=TPM_B_plasma$ensembl_gene_id ,"TPM"=TPM_B_plasma$TPM)
# TPM_B_plasma <- ddply(TPM_B_plasma,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_B_plasma <- merge(TPM_B_plasma,TPM_B_plasma_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_B_plasma$TPM <- log10(TPM_B_plasma$TPM)
# 
# TPM_B_plasma[is.na(TPM_B_plasma)]= 0
# TPM_B_plasma <- subset(TPM_B_plasma, TPM_B_plasma$TPM>=-1)
# print(dim(TPM_B_plasma)) # 13729  2682
# 
# rownames(TPM_B_plasma) <- TPM_B_plasma$ID
# TPM_B_plasma$ID <- NULL
# TPM_B_plasma$ensembl_gene_id <- NULL
# TPM_B_plasma$gene_biotype <- NULL
# TPM_B_plasma$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_B_plasma[is.na(TPM_B_plasma)]=0
# TPM_B_plasma <- TPM_B_plasma[, -caret::nearZeroVar(TPM_B_plasma, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_B_plasma)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_B_plasma), 0.8 * nrow(TPM_B_plasma))
# TPM_B_plasma_train <- TPM_B_plasma[train_row,]
# TPM_B_plasma_test <- TPM_B_plasma[-train_row,]
# 
# dim(TPM_B_plasma_train)
# dim(TPM_B_plasma_test)
# 
# write.table(TPM_B_plasma_train,"~/Analysis/split_dataset/TPM_B_plasma_train",sep=";")
# write.table(TPM_B_plasma_test,"~/Analysis/split_dataset/TPM_B_plasma_test",sep=";")
# 
# TPM_B_plasma_train <- read.delim("~/Analysis/split_dataset/TPM_B_plasma_train",sep=";")
# TPM_B_plasma_test <- read.delim("~/Analysis/split_dataset/TPM_B_plasma_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_B_plasma <- train(TPM~.,
#                         data=TPM_B_plasma_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_B_plasma)
# 
# saveRDS(xgb_TPM_B_plasma, "~/Analysis/models/xgb_TPM_B_plasma_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_B_plasma_03_05_2022.RDS")
# 
# xgb_TPM_B_plasma_predict_test <- predict(xgb_TPM_B_plasma,TPM_B_plasma_test)
# xgb_TPM_B_plasma_predict_test <- data.frame(xgb_TPM_B_plasma_predict_test)
# print(lm(xgb_TPM_B_plasma_predict_test$xgb_TPM_B_plasma_predict_test~TPM_B_plasma_test$TPM))
# print(cor(xgb_TPM_B_plasma_predict_test$xgb_TPM_B_plasma_predict_test,TPM_B_plasma_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_B_plasma,TPM_B_plasma_param)










###__________________________________________________________________________________________###
###_________________________________________ B naive ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for B_S_mem")
# TPM_B_S_mem <- read.delim("~/Analysis/data/immune_cells_RNA/B_S_mem_TPM_log10.csv",sep = ";",dec = ",")
# TPM_B_S_mem <- merge(TPM_B_S_mem,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_B_S_mem_param <- ddply(TPM_B_S_mem,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_B_S_mem_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_B_S_mem$TPM <- 10^(TPM_B_S_mem$TPM)
# TPM_B_S_mem <- data.frame("ID"=TPM_B_S_mem$ensembl_gene_id ,"TPM"=TPM_B_S_mem$TPM)
# TPM_B_S_mem <- ddply(TPM_B_S_mem,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_B_S_mem <- merge(TPM_B_S_mem,TPM_B_S_mem_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_B_S_mem$TPM <- log10(TPM_B_S_mem$TPM)
# 
# TPM_B_S_mem[is.na(TPM_B_S_mem)]= 0
# TPM_B_S_mem <- subset(TPM_B_S_mem, TPM_B_S_mem$TPM>=-1)
# print(dim(TPM_B_S_mem)) # 13729  2682
# 
# rownames(TPM_B_S_mem) <- TPM_B_S_mem$ID
# TPM_B_S_mem$ID <- NULL
# TPM_B_S_mem$ensembl_gene_id <- NULL
# TPM_B_S_mem$gene_biotype <- NULL
# TPM_B_S_mem$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_B_S_mem[is.na(TPM_B_S_mem)]=0
# TPM_B_S_mem <- TPM_B_S_mem[, -caret::nearZeroVar(TPM_B_S_mem, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_B_S_mem)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_B_S_mem), 0.8 * nrow(TPM_B_S_mem))
# TPM_B_S_mem_train <- TPM_B_S_mem[train_row,]
# TPM_B_S_mem_test <- TPM_B_S_mem[-train_row,]
# 
# dim(TPM_B_S_mem_train)
# dim(TPM_B_S_mem_test)
# 
# write.table(TPM_B_S_mem_train,"~/Analysis/split_dataset/TPM_B_S_mem_train",sep=";")
# write.table(TPM_B_S_mem_test,"~/Analysis/split_dataset/TPM_B_S_mem_test",sep=";")
# 
# # TPM_B_S_mem_train <- read.delim("~/Analysis/split_dataset/TPM_B_S_mem_train",sep=";")
# # TPM_B_S_mem_test <- read.delim("~/Analysis/split_dataset/TPM_B_S_mem_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_B_S_mem <- train(TPM~.,
#                         data=TPM_B_S_mem_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_B_S_mem)
# 
# saveRDS(xgb_TPM_B_S_mem, "~/Analysis/models/xgb_TPM_B_S_mem_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_B_S_mem_03_05_2022.RDS")
# 
# xgb_TPM_B_S_mem_predict_test <- predict(xgb_TPM_B_S_mem,TPM_B_S_mem_test)
# xgb_TPM_B_S_mem_predict_test <- data.frame(xgb_TPM_B_S_mem_predict_test)
# print(lm(xgb_TPM_B_S_mem_predict_test$xgb_TPM_B_S_mem_predict_test~TPM_B_S_mem_test$TPM))
# print(cor(xgb_TPM_B_S_mem_predict_test$xgb_TPM_B_S_mem_predict_test,TPM_B_S_mem_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_B_S_mem,TPM_B_S_mem_param)













###__________________________________________________________________________________________###
###_________________________________________ Basophils ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Basophils")
# TPM_Basophils <- read.delim("~/Analysis/data/immune_cells_RNA/Basophils_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Basophils <- merge(TPM_Basophils,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Basophils_param <- ddply(TPM_Basophils,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Basophils_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Basophils$TPM <- 10^(TPM_Basophils$TPM)
# TPM_Basophils <- data.frame("ID"=TPM_Basophils$ensembl_gene_id ,"TPM"=TPM_Basophils$TPM)
# TPM_Basophils <- ddply(TPM_Basophils,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Basophils <- merge(TPM_Basophils,TPM_Basophils_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Basophils$TPM <- log10(TPM_Basophils$TPM)
# 
# TPM_Basophils[is.na(TPM_Basophils)]= 0
# TPM_Basophils <- subset(TPM_Basophils, TPM_Basophils$TPM>=-1)
# print(dim(TPM_Basophils)) # 13729  2682
# 
# rownames(TPM_Basophils) <- TPM_Basophils$ID
# TPM_Basophils$ID <- NULL
# TPM_Basophils$ensembl_gene_id <- NULL
# TPM_Basophils$gene_biotype <- NULL
# TPM_Basophils$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Basophils[is.na(TPM_Basophils)]=0
# TPM_Basophils <- TPM_Basophils[, -caret::nearZeroVar(TPM_Basophils, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Basophils)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Basophils), 0.8 * nrow(TPM_Basophils))
# TPM_Basophils_train <- TPM_Basophils[train_row,]
# TPM_Basophils_test <- TPM_Basophils[-train_row,]
# 
# dim(TPM_Basophils_train)
# dim(TPM_Basophils_test)
# 
# write.table(TPM_Basophils_train,"~/Analysis/split_dataset/TPM_Basophils_train",sep=";")
# write.table(TPM_Basophils_test,"~/Analysis/split_dataset/TPM_Basophils_test",sep=";")
# 
# # TPM_Basophils_train <- read.delim("~/Analysis/split_dataset/TPM_Basophils_train",sep=";")
# # TPM_Basophils_test <- read.delim("~/Analysis/split_dataset/TPM_Basophils_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Basophils <- train(TPM~.,
#                         data=TPM_Basophils_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Basophils)
# 
# saveRDS(xgb_TPM_Basophils, "~/Analysis/models/xgb_TPM_Basophils_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_Basophils_03_05_2022.RDS")
# 
# xgb_TPM_Basophils_predict_test <- predict(xgb_TPM_Basophils,TPM_Basophils_test)
# xgb_TPM_Basophils_predict_test <- data.frame(xgb_TPM_Basophils_predict_test)
# print(lm(xgb_TPM_Basophils_predict_test$xgb_TPM_Basophils_predict_test~TPM_Basophils_test$TPM))
# print(cor(xgb_TPM_Basophils_predict_test$xgb_TPM_Basophils_predict_test,TPM_Basophils_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_Basophils,TPM_Basophils_param)













###__________________________________________________________________________________________###
###_________________________________________ mDC ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for mDC")
# TPM_mDC <- read.delim("~/Analysis/data/immune_cells_RNA/mDC_TPM_log10.csv",sep = ";",dec = ",")
# TPM_mDC <- merge(TPM_mDC,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_mDC_param <- ddply(TPM_mDC,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_mDC_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_mDC$TPM <- 10^(TPM_mDC$TPM)
# TPM_mDC <- data.frame("ID"=TPM_mDC$ensembl_gene_id ,"TPM"=TPM_mDC$TPM)
# TPM_mDC <- ddply(TPM_mDC,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_mDC <- merge(TPM_mDC,TPM_mDC_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_mDC$TPM <- log10(TPM_mDC$TPM)
# 
# TPM_mDC[is.na(TPM_mDC)]= 0
# TPM_mDC <- subset(TPM_mDC, TPM_mDC$TPM>=-1)
# print(dim(TPM_mDC)) # 13729  2682
# 
# rownames(TPM_mDC) <- TPM_mDC$ID
# TPM_mDC$ID <- NULL
# TPM_mDC$ensembl_gene_id <- NULL
# TPM_mDC$gene_biotype <- NULL
# TPM_mDC$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_mDC[is.na(TPM_mDC)]=0
# TPM_mDC <- TPM_mDC[, -caret::nearZeroVar(TPM_mDC, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_mDC)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_mDC), 0.8 * nrow(TPM_mDC))
# TPM_mDC_train <- TPM_mDC[train_row,]
# TPM_mDC_test <- TPM_mDC[-train_row,]
# 
# dim(TPM_mDC_train)
# dim(TPM_mDC_test)
# 
# write.table(TPM_mDC_train,"~/Analysis/split_dataset/TPM_mDC_train",sep=";")
# write.table(TPM_mDC_test,"~/Analysis/split_dataset/TPM_mDC_test",sep=";")
# 
# # TPM_mDC_train <- read.delim("~/Analysis/split_dataset/TPM_mDC_train",sep=";")
# # TPM_mDC_test <- read.delim("~/Analysis/split_dataset/TPM_mDC_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_mDC <- train(TPM~.,
#                         data=TPM_mDC_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_mDC)
# 
# saveRDS(xgb_TPM_mDC, "~/Analysis/models/xgb_TPM_mDC_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_mDC_03_05_2022.RDS")
# 
# xgb_TPM_mDC_predict_test <- predict(xgb_TPM_mDC,TPM_mDC_test)
# xgb_TPM_mDC_predict_test <- data.frame(xgb_TPM_mDC_predict_test)
# print(lm(xgb_TPM_mDC_predict_test$xgb_TPM_mDC_predict_test~TPM_mDC_test$TPM))
# print(cor(xgb_TPM_mDC_predict_test$xgb_TPM_mDC_predict_test,TPM_mDC_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_mDC,TPM_mDC_param)













###__________________________________________________________________________________________###
###_________________________________________ Mono_classical ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Mono_classical")
# TPM_Mono_classical <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_classical_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Mono_classical <- merge(TPM_Mono_classical,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Mono_classical_param <- ddply(TPM_Mono_classical,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Mono_classical_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Mono_classical$TPM <- 10^(TPM_Mono_classical$TPM)
# TPM_Mono_classical <- data.frame("ID"=TPM_Mono_classical$ensembl_gene_id ,"TPM"=TPM_Mono_classical$TPM)
# TPM_Mono_classical <- ddply(TPM_Mono_classical,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Mono_classical <- merge(TPM_Mono_classical,TPM_Mono_classical_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Mono_classical$TPM <- log10(TPM_Mono_classical$TPM)
# 
# TPM_Mono_classical[is.na(TPM_Mono_classical)]= 0
# TPM_Mono_classical <- subset(TPM_Mono_classical, TPM_Mono_classical$TPM>=-1)
# print(dim(TPM_Mono_classical)) # 13729  2682
# 
# rownames(TPM_Mono_classical) <- TPM_Mono_classical$ID
# TPM_Mono_classical$ID <- NULL
# TPM_Mono_classical$ensembl_gene_id <- NULL
# TPM_Mono_classical$gene_biotype <- NULL
# TPM_Mono_classical$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Mono_classical[is.na(TPM_Mono_classical)]=0
# TPM_Mono_classical <- TPM_Mono_classical[, -caret::nearZeroVar(TPM_Mono_classical, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Mono_classical)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Mono_classical), 0.8 * nrow(TPM_Mono_classical))
# TPM_Mono_classical_train <- TPM_Mono_classical[train_row,]
# TPM_Mono_classical_test <- TPM_Mono_classical[-train_row,]
# 
# dim(TPM_Mono_classical_train)
# dim(TPM_Mono_classical_test)
# 
# write.table(TPM_Mono_classical_train,"~/Analysis/split_dataset/TPM_Mono_classical_train",sep=";")
# write.table(TPM_Mono_classical_test,"~/Analysis/split_dataset/TPM_Mono_classical_test",sep=";")
# 
# # TPM_Mono_classical_train <- read.delim("~/Analysis/split_dataset/TPM_Mono_classical_train",sep=";")
# # TPM_Mono_classical_test <- read.delim("~/Analysis/split_dataset/TPM_Mono_classical_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Mono_classical <- train(TPM~.,
#                         data=TPM_Mono_classical_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Mono_classical)
# 
# saveRDS(xgb_TPM_Mono_classical, "~/Analysis/models/xgb_TPM_Mono_classical_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_Mono_classical_03_05_2022.RDS")
# 
# xgb_TPM_Mono_classical_predict_test <- predict(xgb_TPM_Mono_classical,TPM_Mono_classical_test)
# xgb_TPM_Mono_classical_predict_test <- data.frame(xgb_TPM_Mono_classical_predict_test)
# print(lm(xgb_TPM_Mono_classical_predict_test$xgb_TPM_Mono_classical_predict_test~TPM_Mono_classical_test$TPM))
# print(cor(xgb_TPM_Mono_classical_predict_test$xgb_TPM_Mono_classical_predict_test,TPM_Mono_classical_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_Mono_classical,TPM_Mono_classical_param)













###__________________________________________________________________________________________###
###_________________________________________ Mono_inter ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Mono_inter")
# TPM_Mono_inter <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_inter_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Mono_inter <- merge(TPM_Mono_inter,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Mono_inter_param <- ddply(TPM_Mono_inter,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Mono_inter_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Mono_inter$TPM <- 10^(TPM_Mono_inter$TPM)
# TPM_Mono_inter <- data.frame("ID"=TPM_Mono_inter$ensembl_gene_id ,"TPM"=TPM_Mono_inter$TPM)
# TPM_Mono_inter <- ddply(TPM_Mono_inter,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Mono_inter <- merge(TPM_Mono_inter,TPM_Mono_inter_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Mono_inter$TPM <- log10(TPM_Mono_inter$TPM)
# 
# TPM_Mono_inter[is.na(TPM_Mono_inter)]= 0
# TPM_Mono_inter <- subset(TPM_Mono_inter, TPM_Mono_inter$TPM>=-1)
# print(dim(TPM_Mono_inter)) # 13729  2682
# 
# rownames(TPM_Mono_inter) <- TPM_Mono_inter$ID
# TPM_Mono_inter$ID <- NULL
# TPM_Mono_inter$ensembl_gene_id <- NULL
# TPM_Mono_inter$gene_biotype <- NULL
# TPM_Mono_inter$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Mono_inter[is.na(TPM_Mono_inter)]=0
# TPM_Mono_inter <- TPM_Mono_inter[, -caret::nearZeroVar(TPM_Mono_inter, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Mono_inter)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Mono_inter), 0.8 * nrow(TPM_Mono_inter))
# TPM_Mono_inter_train <- TPM_Mono_inter[train_row,]
# TPM_Mono_inter_test <- TPM_Mono_inter[-train_row,]
# 
# dim(TPM_Mono_inter_train)
# dim(TPM_Mono_inter_test)
# 
# write.table(TPM_Mono_inter_train,"~/Analysis/split_dataset/TPM_Mono_inter_train",sep=";")
# write.table(TPM_Mono_inter_test,"~/Analysis/split_dataset/TPM_Mono_inter_test",sep=";")
# 
# # TPM_Mono_inter_train <- read.delim("~/Analysis/split_dataset/TPM_Mono_inter_train",sep=";")
# # TPM_Mono_inter_test <- read.delim("~/Analysis/split_dataset/TPM_Mono_inter_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000, 
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Mono_inter <- train(TPM~.,
#                         data=TPM_Mono_inter_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Mono_inter)
# 
# saveRDS(xgb_TPM_Mono_inter, "~/Analysis/models/xgb_TPM_Mono_inter_03_05_2022.RDS")
# 
# # xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_Mono_inter_03_05_2022.RDS")
# 
# xgb_TPM_Mono_inter_predict_test <- predict(xgb_TPM_Mono_inter,TPM_Mono_inter_test)
# xgb_TPM_Mono_inter_predict_test <- data.frame(xgb_TPM_Mono_inter_predict_test)
# print(lm(xgb_TPM_Mono_inter_predict_test$xgb_TPM_Mono_inter_predict_test~TPM_Mono_inter_test$TPM))
# print(cor(xgb_TPM_Mono_inter_predict_test$xgb_TPM_Mono_inter_predict_test,TPM_Mono_inter_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_Mono_inter,TPM_Mono_inter_param)









### RESTART HERE 06-05-22 ###












###__________________________________________________________________________________________###
###_________________________________________ Mono_non_classical ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Mono_non_classical")
# TPM_Mono_non_classical <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_non_classical_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Mono_non_classical <- merge(TPM_Mono_non_classical,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Mono_non_classical_param <- ddply(TPM_Mono_non_classical,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Mono_non_classical_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Mono_non_classical$TPM <- 10^(TPM_Mono_non_classical$TPM)
# TPM_Mono_non_classical <- data.frame("ID"=TPM_Mono_non_classical$ensembl_gene_id ,"TPM"=TPM_Mono_non_classical$TPM)
# TPM_Mono_non_classical <- ddply(TPM_Mono_non_classical,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Mono_non_classical <- merge(TPM_Mono_non_classical,TPM_Mono_non_classical_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Mono_non_classical$TPM <- log10(TPM_Mono_non_classical$TPM)
# 
# TPM_Mono_non_classical[is.na(TPM_Mono_non_classical)]= 0
# TPM_Mono_non_classical <- subset(TPM_Mono_non_classical, TPM_Mono_non_classical$TPM>=-1)
# print(dim(TPM_Mono_non_classical)) # 13729  2682
# 
# rownames(TPM_Mono_non_classical) <- TPM_Mono_non_classical$ID
# TPM_Mono_non_classical$ID <- NULL
# TPM_Mono_non_classical$ensembl_gene_id <- NULL
# TPM_Mono_non_classical$gene_biotype <- NULL
# TPM_Mono_non_classical$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Mono_non_classical[is.na(TPM_Mono_non_classical)]=0
# TPM_Mono_non_classical <- TPM_Mono_non_classical[, -caret::nearZeroVar(TPM_Mono_non_classical, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Mono_non_classical)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Mono_non_classical), 0.8 * nrow(TPM_Mono_non_classical))
# TPM_Mono_non_classical_train <- TPM_Mono_non_classical[train_row,]
# TPM_Mono_non_classical_test <- TPM_Mono_non_classical[-train_row,]
# 
# dim(TPM_Mono_non_classical_train)
# dim(TPM_Mono_non_classical_test)
# 
# write.table(TPM_Mono_non_classical_train,"~/Analysis/split_dataset/TPM_Mono_non_classical_train",sep=";")
# write.table(TPM_Mono_non_classical_test,"~/Analysis/split_dataset/TPM_Mono_non_classical_test",sep=";")
# 
# TPM_Mono_non_classical_train <- read.delim("~/Analysis/split_dataset/TPM_Mono_non_classical_train",sep=";")
# TPM_Mono_non_classical_test <- read.delim("~/Analysis/split_dataset/TPM_Mono_non_classical_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000,
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Mono_non_classical <- train(TPM~.,
#                         data=TPM_Mono_non_classical_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread= 60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Mono_non_classical)
# 
# saveRDS(xgb_TPM_Mono_non_classical, "~/Analysis/models/xgb_TPM_Mono_non_classical_03_05_2022.RDS")
# 
# xgb_TPM_Mono_non_classical <- readRDS( "~/Analysis/models/xgb_TPM_Mono_non_classical_03_05_2022.RDS")
# # 
# xgb_TPM_Mono_non_classical_predict_test <- predict(xgb_TPM_Mono_non_classical,TPM_Mono_non_classical_test)
# xgb_TPM_Mono_non_classical_predict_test <- data.frame(xgb_TPM_Mono_non_classical_predict_test)
# print(lm(xgb_TPM_Mono_non_classical_predict_test$xgb_TPM_Mono_non_classical_predict_test~TPM_Mono_non_classical_test$TPM))
# print(cor(xgb_TPM_Mono_non_classical_predict_test$xgb_TPM_Mono_non_classical_predict_test,TPM_Mono_non_classical_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_Mono_non_classical,TPM_Mono_non_classical_param)

















###__________________________________________________________________________________________###
###_________________________________________ Neutrophils ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Neutrophils")
# TPM_Neutrophils <- read.delim("~/Analysis/data/immune_cells_RNA/Neutrophils_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Neutrophils <- merge(TPM_Neutrophils,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Neutrophils_param <- ddply(TPM_Neutrophils,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Neutrophils_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Neutrophils$TPM <- 10^(TPM_Neutrophils$TPM)
# TPM_Neutrophils <- data.frame("ID"=TPM_Neutrophils$ensembl_gene_id ,"TPM"=TPM_Neutrophils$TPM)
# TPM_Neutrophils <- ddply(TPM_Neutrophils,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Neutrophils <- merge(TPM_Neutrophils,TPM_Neutrophils_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Neutrophils$TPM <- log10(TPM_Neutrophils$TPM)
# 
# TPM_Neutrophils[is.na(TPM_Neutrophils)]= 0
# TPM_Neutrophils <- subset(TPM_Neutrophils, TPM_Neutrophils$TPM>=-1)
# print(dim(TPM_Neutrophils)) # 13729  2682
# 
# rownames(TPM_Neutrophils) <- TPM_Neutrophils$ID
# TPM_Neutrophils$ID <- NULL
# TPM_Neutrophils$ensembl_gene_id <- NULL
# TPM_Neutrophils$gene_biotype <- NULL
# TPM_Neutrophils$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Neutrophils[is.na(TPM_Neutrophils)]=0
# TPM_Neutrophils <- TPM_Neutrophils[, -caret::nearZeroVar(TPM_Neutrophils, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Neutrophils)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Neutrophils), 0.8 * nrow(TPM_Neutrophils))
# TPM_Neutrophils_train <- TPM_Neutrophils[train_row,]
# TPM_Neutrophils_test <- TPM_Neutrophils[-train_row,]
# 
# dim(TPM_Neutrophils_train)
# dim(TPM_Neutrophils_test)
# 
# write.table(TPM_Neutrophils_train,"~/Analysis/split_dataset/TPM_Neutrophils_train",sep=";")
# write.table(TPM_Neutrophils_test,"~/Analysis/split_dataset/TPM_Neutrophils_test",sep=";")
# 
# TPM_Neutrophils_train <- read.delim("~/Analysis/split_dataset/TPM_Neutrophils_train",sep=";")
# TPM_Neutrophils_test <- read.delim("~/Analysis/split_dataset/TPM_Neutrophils_test",sep=";")

# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000,
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Neutrophils <- train(TPM~.,
#                         data=TPM_Neutrophils_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Neutrophils)
# 
# saveRDS(xgb_TPM_Neutrophils, "~/Analysis/models/xgb_TPM_Neutrophils_03_05_2022.RDS")
# 
# xgb_TPM_Neutrophils <- readRDS( "~/Analysis/models/xgb_TPM_Neutrophils_03_05_2022.RDS")
# 
# xgb_TPM_Neutrophils_predict_test <- predict(xgb_TPM_Neutrophils,TPM_Neutrophils_test)
# xgb_TPM_Neutrophils_predict_test <- data.frame(xgb_TPM_Neutrophils_predict_test)
# print(lm(xgb_TPM_Neutrophils_predict_test$xgb_TPM_Neutrophils_predict_test~TPM_Neutrophils_test$TPM))
# print(cor(xgb_TPM_Neutrophils_predict_test$xgb_TPM_Neutrophils_predict_test,TPM_Neutrophils_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
plot(xgb_TPM_Neutrophils_predict_test$xgb_TPM_Neutrophils_predict_test,TPM_Neutrophils_test$TPM)


# gc()
# rm(TPM_Neutrophils,TPM_Neutrophils_param)














###__________________________________________________________________________________________###
###_________________________________________ NK ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for NK")
# TPM_NK <- read.delim("~/Analysis/data/immune_cells_RNA/NK_TPM_log10.csv",sep = ";",dec = ",")
# TPM_NK <- merge(TPM_NK,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_NK_param <- ddply(TPM_NK,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_NK_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_NK$TPM <- 10^(TPM_NK$TPM)
# TPM_NK <- data.frame("ID"=TPM_NK$ensembl_gene_id ,"TPM"=TPM_NK$TPM)
# TPM_NK <- ddply(TPM_NK,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_NK <- merge(TPM_NK,TPM_NK_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_NK$TPM <- log10(TPM_NK$TPM)
# 
# TPM_NK[is.na(TPM_NK)]= 0
# TPM_NK <- subset(TPM_NK, TPM_NK$TPM>=-1)
# print(dim(TPM_NK)) # 13729  2682
# 
# rownames(TPM_NK) <- TPM_NK$ID
# TPM_NK$ID <- NULL
# TPM_NK$ensembl_gene_id <- NULL
# TPM_NK$gene_biotype <- NULL
# TPM_NK$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_NK[is.na(TPM_NK)]=0
# TPM_NK <- TPM_NK[, -caret::nearZeroVar(TPM_NK, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_NK)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_NK), 0.8 * nrow(TPM_NK))
# TPM_NK_train <- TPM_NK[train_row,]
# TPM_NK_test <- TPM_NK[-train_row,]
# 
# dim(TPM_NK_train)
# dim(TPM_NK_test)
# 
# write.table(TPM_NK_train,"~/Analysis/split_dataset/TPM_NK_train",sep=";")
# write.table(TPM_NK_test,"~/Analysis/split_dataset/TPM_NK_test",sep=";")
# 
# TPM_NK_train <- read.delim("~/Analysis/split_dataset/TPM_NK_train",sep=";")
# TPM_NK_test <- read.delim("~/Analysis/split_dataset/TPM_NK_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000,
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_NK <- train(TPM~.,
#                         data=TPM_NK_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_NK)
# 
# saveRDS(xgb_TPM_NK, "~/Analysis/models/xgb_TPM_NK_03_05_2022.RDS")
# 
# xgb_TPM_NK <- readRDS( "~/Analysis/models/xgb_TPM_NK_03_05_2022.RDS")
# 
# xgb_TPM_NK_predict_test <- predict(xgb_TPM_NK,TPM_NK_test)
# xgb_TPM_NK_predict_test <- data.frame(xgb_TPM_NK_predict_test)
# print(lm(xgb_TPM_NK_predict_test$xgb_TPM_NK_predict_test~TPM_NK_test$TPM))
# print(cor(xgb_TPM_NK_predict_test$xgb_TPM_NK_predict_test,TPM_NK_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_NK,TPM_NK_param)














###__________________________________________________________________________________________###
###_________________________________________ pDC ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for pDC")
# TPM_pDC <- read.delim("~/Analysis/data/immune_cells_RNA/pDC_TPM_log10.csv",sep = ";",dec = ",")
# TPM_pDC <- merge(TPM_pDC,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_pDC_param <- ddply(TPM_pDC,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_pDC_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_pDC$TPM <- 10^(TPM_pDC$TPM)
# TPM_pDC <- data.frame("ID"=TPM_pDC$ensembl_gene_id ,"TPM"=TPM_pDC$TPM)
# TPM_pDC <- ddply(TPM_pDC,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_pDC <- merge(TPM_pDC,TPM_pDC_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_pDC$TPM <- log10(TPM_pDC$TPM)
# 
# TPM_pDC[is.na(TPM_pDC)]= 0
# TPM_pDC <- subset(TPM_pDC, TPM_pDC$TPM>=-1)
# print(dim(TPM_pDC)) # 13729  2682
# 
# rownames(TPM_pDC) <- TPM_pDC$ID
# TPM_pDC$ID <- NULL
# TPM_pDC$ensembl_gene_id <- NULL
# TPM_pDC$gene_biotype <- NULL
# TPM_pDC$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_pDC[is.na(TPM_pDC)]=0
# TPM_pDC <- TPM_pDC[, -caret::nearZeroVar(TPM_pDC, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_pDC)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_pDC), 0.8 * nrow(TPM_pDC))
# TPM_pDC_train <- TPM_pDC[train_row,]
# TPM_pDC_test <- TPM_pDC[-train_row,]
# 
# dim(TPM_pDC_train)
# dim(TPM_pDC_test)
# 
# write.table(TPM_pDC_train,"~/Analysis/split_dataset/TPM_pDC_train",sep=";")
# write.table(TPM_pDC_test,"~/Analysis/split_dataset/TPM_pDC_test",sep=";")
# 
# TPM_pDC_train <- read.delim("~/Analysis/split_dataset/TPM_pDC_train",sep=";")
# TPM_pDC_test <- read.delim("~/Analysis/split_dataset/TPM_pDC_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000,
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_pDC <- train(TPM~.,
#                         data=TPM_pDC_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_pDC)
# 
# saveRDS(xgb_TPM_pDC, "~/Analysis/models/xgb_TPM_pDC_03_05_2022.RDS")
# 
# xgb_TPM_pDC <- readRDS( "~/Analysis/models/xgb_TPM_pDC_03_05_2022.RDS")
# 
# xgb_TPM_pDC_predict_test <- predict(xgb_TPM_pDC,TPM_pDC_test)
# xgb_TPM_pDC_predict_test <- data.frame(xgb_TPM_pDC_predict_test)
# print(lm(xgb_TPM_pDC_predict_test$xgb_TPM_pDC_predict_test~TPM_pDC_test$TPM))
# print(cor(xgb_TPM_pDC_predict_test$xgb_TPM_pDC_predict_test,TPM_pDC_test$TPM, method = "pearson", use = "complete.obs")^2)
# 
# 
# gc()
# rm(TPM_pDC,TPM_pDC_param)

















###__________________________________________________________________________________________###
###_________________________________________ Tregs ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for Tregs")
# TPM_Tregs <- read.delim("~/Analysis/data/immune_cells_RNA/Tregs_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Tregs <- merge(TPM_Tregs,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Tregs_param <- ddply(TPM_Tregs,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Tregs_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Tregs$TPM <- 10^(TPM_Tregs$TPM)
# TPM_Tregs <- data.frame("ID"=TPM_Tregs$ensembl_gene_id ,"TPM"=TPM_Tregs$TPM)
# TPM_Tregs <- ddply(TPM_Tregs,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Tregs <- merge(TPM_Tregs,TPM_Tregs_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Tregs$TPM <- log10(TPM_Tregs$TPM)
# 
# TPM_Tregs[is.na(TPM_Tregs)]= 0
# TPM_Tregs <- subset(TPM_Tregs, TPM_Tregs$TPM>=-1)
# print(dim(TPM_Tregs)) # 13729  2682
# 
# rownames(TPM_Tregs) <- TPM_Tregs$ID
# TPM_Tregs$ID <- NULL
# TPM_Tregs$ensembl_gene_id <- NULL
# TPM_Tregs$gene_biotype <- NULL
# TPM_Tregs$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Tregs[is.na(TPM_Tregs)]=0
# TPM_Tregs <- TPM_Tregs[, -caret::nearZeroVar(TPM_Tregs, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Tregs)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Tregs), 0.8 * nrow(TPM_Tregs))
# TPM_Tregs_train <- TPM_Tregs[train_row,]
# TPM_Tregs_test <- TPM_Tregs[-train_row,]
# 
# dim(TPM_Tregs_train)
# dim(TPM_Tregs_test)
# 
# write.table(TPM_Tregs_train,"~/Analysis/split_dataset/TPM_Tregs_train",sep=";")
# write.table(TPM_Tregs_test,"~/Analysis/split_dataset/TPM_Tregs_test",sep=";")
# 
# TPM_Tregs_train <- read.delim("~/Analysis/split_dataset/TPM_Tregs_train",sep=";")
# TPM_Tregs_test <- read.delim("~/Analysis/split_dataset/TPM_Tregs_test",sep=";")
# 
# 
# 
# ## Training the XGB model ##
# print("modeling")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid <- expand.grid(nrounds = 10000,
#                        max_depth = 6,
#                        colsample_bytree = 0.4,
#                        eta = 0.05,
#                        gamma=0,
#                        min_child_weight = 0.9,
#                        subsample = 1)
# 
# start_time <- Sys.time()
# xgb_TPM_Tregs <- train(TPM~.,
#                         data=TPM_Tregs_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid,
#                         na.action = na.omit,
#                         nthread=60,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of  hours
# 
# print(xgb_TPM_Tregs)
# 
# saveRDS(xgb_TPM_Tregs, "~/Analysis/models/xgb_TPM_Tregs_03_05_2022.RDS")
# 
xgb_TPM_Tregs <- readRDS( "~/Analysis/models/xgb_TPM_Tregs_03_05_2022.RDS")

xgb_TPM_Tregs_predict_test <- predict(xgb_TPM_Tregs,TPM_Tregs_test)
xgb_TPM_Tregs_predict_test <- data.frame(xgb_TPM_Tregs_predict_test)
print(lm(xgb_TPM_Tregs_predict_test$xgb_TPM_Tregs_predict_test~TPM_Tregs_test$TPM))
print(cor(xgb_TPM_Tregs_predict_test$xgb_TPM_Tregs_predict_test,TPM_Tregs_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_Tregs_predict_test$xgb_TPM_Tregs_predict_test,TPM_Tregs_test$TPM)

# 
# 
# gc()
# rm(TPM_Tregs,TPM_Tregs_param)













































###__________________________________________________________________________________________###
###_________________________________________ NT ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for NT")
TPM_NT <- read.delim("~/Analysis/RBP_validation/data/NT_TPM_log10.csv",sep = ";",dec = ",")
TPM_NT <- merge(TPM_NT,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_NT_param <- ddply(TPM_NT,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_NT_param$TPM <- NULL

print("integrating TPM per gene")
TPM_NT$TPM <- 10^(TPM_NT$TPM)
TPM_NT <- data.frame("ID"=TPM_NT$ensembl_gene_id ,"TPM"=TPM_NT$TPM)
TPM_NT <- ddply(TPM_NT,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_NT <- merge(TPM_NT,TPM_NT_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_NT$TPM <- log10(TPM_NT$TPM)

TPM_NT[is.na(TPM_NT)]= 0
TPM_NT <- subset(TPM_NT, TPM_NT$TPM>=-1)
print(dim(TPM_NT)) # 13729  2682

rownames(TPM_NT) <- TPM_NT$ID
TPM_NT$ID <- NULL
TPM_NT$ensembl_gene_id <- NULL
TPM_NT$gene_biotype <- NULL
TPM_NT$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_NT[is.na(TPM_NT)]=0
TPM_NT <- TPM_NT[, -caret::nearZeroVar(TPM_NT, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_NT)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_NT), 0.8 * nrow(TPM_NT))
TPM_NT_train <- TPM_NT[train_row,]
TPM_NT_test <- TPM_NT[-train_row,]

dim(TPM_NT_train)
dim(TPM_NT_test)

write.table(TPM_NT_train,"~/Analysis/RBP_validation/split_datasets/TPM_NT_train",sep=";")
write.table(TPM_NT_test,"~/Analysis/RBP_validation/split_datasets/TPM_NT_test",sep=";")

# TPM_NT_train <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_NT_train",sep=";")
# TPM_NT_test <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_NT_test",sep=";")



## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)

start_time <- Sys.time()
xgb_TPM_NT <- train(TPM~.,
                        data=TPM_NT_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=60,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_NT)

saveRDS(xgb_TPM_NT, "~/Analysis/RBP_validation/models/xgb_TPM_NT_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/RBP_validation/models/xgb_TPM_NT_03_05_2022.RDS")

xgb_TPM_NT_predict_test <- predict(xgb_TPM_NT,TPM_NT_test)
xgb_TPM_NT_predict_test <- data.frame(xgb_TPM_NT_predict_test)
print(lm(xgb_TPM_NT_predict_test$xgb_TPM_NT_predict_test~TPM_NT_test$TPM))
print(cor(xgb_TPM_NT_predict_test$xgb_TPM_NT_predict_test,TPM_NT_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_NT,TPM_NT_param)



















###__________________________________________________________________________________________###
###_________________________________________ ELAVL1 ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for ELAVL1")
TPM_ELAVL1 <- read.delim("~/Analysis/RBP_validation/data/ELAVL1_TPM_log10.csv",sep = ";",dec = ",")
TPM_ELAVL1 <- merge(TPM_ELAVL1,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_ELAVL1_param <- ddply(TPM_ELAVL1,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_ELAVL1_param$TPM <- NULL

print("integrating TPM per gene")
TPM_ELAVL1$TPM <- 10^(TPM_ELAVL1$TPM)
TPM_ELAVL1 <- data.frame("ID"=TPM_ELAVL1$ensembl_gene_id ,"TPM"=TPM_ELAVL1$TPM)
TPM_ELAVL1 <- ddply(TPM_ELAVL1,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_ELAVL1 <- merge(TPM_ELAVL1,TPM_ELAVL1_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_ELAVL1$TPM <- log10(TPM_ELAVL1$TPM)

TPM_ELAVL1[is.na(TPM_ELAVL1)]= 0
TPM_ELAVL1 <- subset(TPM_ELAVL1, TPM_ELAVL1$TPM>=-1)
print(dim(TPM_ELAVL1)) # 13729  2682

rownames(TPM_ELAVL1) <- TPM_ELAVL1$ID
TPM_ELAVL1$ID <- NULL
TPM_ELAVL1$ensembl_gene_id <- NULL
TPM_ELAVL1$gene_biotype <- NULL
TPM_ELAVL1$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_ELAVL1[is.na(TPM_ELAVL1)]=0
TPM_ELAVL1 <- TPM_ELAVL1[, -caret::nearZeroVar(TPM_ELAVL1, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_ELAVL1)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_ELAVL1), 0.8 * nrow(TPM_ELAVL1))
TPM_ELAVL1_train <- TPM_ELAVL1[train_row,]
TPM_ELAVL1_test <- TPM_ELAVL1[-train_row,]

dim(TPM_ELAVL1_train)
dim(TPM_ELAVL1_test)

write.table(TPM_ELAVL1_train,"~/Analysis/RBP_validation/split_datasets/TPM_ELAVL1_train",sep=";")
write.table(TPM_ELAVL1_test,"~/Analysis/RBP_validation/split_datasets/TPM_ELAVL1_test",sep=";")

# TPM_ELAVL1_train <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_ELAVL1_train",sep=";")
# TPM_ELAVL1_test <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_ELAVL1_test",sep=";")



## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)

start_time <- Sys.time()
xgb_TPM_ELAVL1 <- train(TPM~.,
                        data=TPM_ELAVL1_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=60,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_ELAVL1)

saveRDS(xgb_TPM_ELAVL1, "~/Analysis/RBP_validation/models/xgb_TPM_ELAVL1_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/RBP_validation/models/xgb_TPM_ELAVL1_03_05_2022.RDS")

xgb_TPM_ELAVL1_predict_test <- predict(xgb_TPM_ELAVL1,TPM_ELAVL1_test)
xgb_TPM_ELAVL1_predict_test <- data.frame(xgb_TPM_ELAVL1_predict_test)
print(lm(xgb_TPM_ELAVL1_predict_test$xgb_TPM_ELAVL1_predict_test~TPM_ELAVL1_test$TPM))
print(cor(xgb_TPM_ELAVL1_predict_test$xgb_TPM_ELAVL1_predict_test,TPM_ELAVL1_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_ELAVL1,TPM_ELAVL1_param)
















###__________________________________________________________________________________________###
###_________________________________________ FUBP3 ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for FUBP3")
TPM_FUBP3 <- read.delim("~/Analysis/RBP_validation/data/FUBP3_TPM_log10.csv",sep = ";",dec = ",")
TPM_FUBP3 <- merge(TPM_FUBP3,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_FUBP3_param <- ddply(TPM_FUBP3,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_FUBP3_param$TPM <- NULL

print("integrating TPM per gene")
TPM_FUBP3$TPM <- 10^(TPM_FUBP3$TPM)
TPM_FUBP3 <- data.frame("ID"=TPM_FUBP3$ensembl_gene_id ,"TPM"=TPM_FUBP3$TPM)
TPM_FUBP3 <- ddply(TPM_FUBP3,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_FUBP3 <- merge(TPM_FUBP3,TPM_FUBP3_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_FUBP3$TPM <- log10(TPM_FUBP3$TPM)

TPM_FUBP3[is.na(TPM_FUBP3)]= 0
TPM_FUBP3 <- subset(TPM_FUBP3, TPM_FUBP3$TPM>=-1)
print(dim(TPM_FUBP3)) # 13729  2682

rownames(TPM_FUBP3) <- TPM_FUBP3$ID
TPM_FUBP3$ID <- NULL
TPM_FUBP3$ensembl_gene_id <- NULL
TPM_FUBP3$gene_biotype <- NULL
TPM_FUBP3$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_FUBP3[is.na(TPM_FUBP3)]=0
TPM_FUBP3 <- TPM_FUBP3[, -caret::nearZeroVar(TPM_FUBP3, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_FUBP3)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_FUBP3), 0.8 * nrow(TPM_FUBP3))
TPM_FUBP3_train <- TPM_FUBP3[train_row,]
TPM_FUBP3_test <- TPM_FUBP3[-train_row,]

dim(TPM_FUBP3_train)
dim(TPM_FUBP3_test)

write.table(TPM_FUBP3_train,"~/Analysis/RBP_validation/split_datasets/TPM_FUBP3_train",sep=";")
write.table(TPM_FUBP3_test,"~/Analysis/RBP_validation/split_datasets/TPM_FUBP3_test",sep=";")

# TPM_FUBP3_train <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_FUBP3_train",sep=";")
# TPM_FUBP3_test <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_FUBP3_test",sep=";")



## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)

start_time <- Sys.time()
xgb_TPM_FUBP3 <- train(TPM~.,
                        data=TPM_FUBP3_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=60,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_FUBP3)

saveRDS(xgb_TPM_FUBP3, "~/Analysis/RBP_validation/models/xgb_TPM_FUBP3_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/RBP_validation/models/xgb_TPM_FUBP3_03_05_2022.RDS")

xgb_TPM_FUBP3_predict_test <- predict(xgb_TPM_FUBP3,TPM_FUBP3_test)
xgb_TPM_FUBP3_predict_test <- data.frame(xgb_TPM_FUBP3_predict_test)
print(lm(xgb_TPM_FUBP3_predict_test$xgb_TPM_FUBP3_predict_test~TPM_FUBP3_test$TPM))
print(cor(xgb_TPM_FUBP3_predict_test$xgb_TPM_FUBP3_predict_test,TPM_FUBP3_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_FUBP3,TPM_FUBP3_param)















###__________________________________________________________________________________________###
###_________________________________________ FXR2 ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for FXR2")
TPM_FXR2 <- read.delim("~/Analysis/RBP_validation/data/FXR2_TPM_log10.csv",sep = ";",dec = ",")
TPM_FXR2 <- merge(TPM_FXR2,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_FXR2_param <- ddply(TPM_FXR2,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_FXR2_param$TPM <- NULL

print("integrating TPM per gene")
TPM_FXR2$TPM <- 10^(TPM_FXR2$TPM)
TPM_FXR2 <- data.frame("ID"=TPM_FXR2$ensembl_gene_id ,"TPM"=TPM_FXR2$TPM)
TPM_FXR2 <- ddply(TPM_FXR2,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_FXR2 <- merge(TPM_FXR2,TPM_FXR2_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_FXR2$TPM <- log10(TPM_FXR2$TPM)

TPM_FXR2[is.na(TPM_FXR2)]= 0
TPM_FXR2 <- subset(TPM_FXR2, TPM_FXR2$TPM>=-1)
print(dim(TPM_FXR2)) # 13729  2682

rownames(TPM_FXR2) <- TPM_FXR2$ID
TPM_FXR2$ID <- NULL
TPM_FXR2$ensembl_gene_id <- NULL
TPM_FXR2$gene_biotype <- NULL
TPM_FXR2$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_FXR2[is.na(TPM_FXR2)]=0
TPM_FXR2 <- TPM_FXR2[, -caret::nearZeroVar(TPM_FXR2, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_FXR2)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_FXR2), 0.8 * nrow(TPM_FXR2))
TPM_FXR2_train <- TPM_FXR2[train_row,]
TPM_FXR2_test <- TPM_FXR2[-train_row,]

dim(TPM_FXR2_train)
dim(TPM_FXR2_test)

write.table(TPM_FXR2_train,"~/Analysis/RBP_validation/split_datasets/TPM_FXR2_train",sep=";")
write.table(TPM_FXR2_test,"~/Analysis/RBP_validation/split_datasets/TPM_FXR2_test",sep=";")

# TPM_FXR2_train <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_FXR2_train",sep=";")
# TPM_FXR2_test <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_FXR2_test",sep=";")



## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)

start_time <- Sys.time()
xgb_TPM_FXR2 <- train(TPM~.,
                        data=TPM_FXR2_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=60,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_FXR2)

saveRDS(xgb_TPM_FXR2, "~/Analysis/RBP_validation/models/xgb_TPM_FXR2_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/RBP_validation/models/xgb_TPM_FXR2_03_05_2022.RDS")

xgb_TPM_FXR2_predict_test <- predict(xgb_TPM_FXR2,TPM_FXR2_test)
xgb_TPM_FXR2_predict_test <- data.frame(xgb_TPM_FXR2_predict_test)
print(lm(xgb_TPM_FXR2_predict_test$xgb_TPM_FXR2_predict_test~TPM_FXR2_test$TPM))
print(cor(xgb_TPM_FXR2_predict_test$xgb_TPM_FXR2_predict_test,TPM_FXR2_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_FXR2,TPM_FXR2_param)















###__________________________________________________________________________________________###
###_________________________________________ PDCD4 ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PDCD4")
TPM_PDCD4 <- read.delim("~/Analysis/RBP_validation/data/PDCD4_TPM_log10.csv",sep = ";",dec = ",")
TPM_PDCD4 <- merge(TPM_PDCD4,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PDCD4_param <- ddply(TPM_PDCD4,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PDCD4_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PDCD4$TPM <- 10^(TPM_PDCD4$TPM)
TPM_PDCD4 <- data.frame("ID"=TPM_PDCD4$ensembl_gene_id ,"TPM"=TPM_PDCD4$TPM)
TPM_PDCD4 <- ddply(TPM_PDCD4,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PDCD4 <- merge(TPM_PDCD4,TPM_PDCD4_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PDCD4$TPM <- log10(TPM_PDCD4$TPM)

TPM_PDCD4[is.na(TPM_PDCD4)]= 0
TPM_PDCD4 <- subset(TPM_PDCD4, TPM_PDCD4$TPM>=-1)
print(dim(TPM_PDCD4)) # 13729  2682

rownames(TPM_PDCD4) <- TPM_PDCD4$ID
TPM_PDCD4$ID <- NULL
TPM_PDCD4$ensembl_gene_id <- NULL
TPM_PDCD4$gene_biotype <- NULL
TPM_PDCD4$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PDCD4[is.na(TPM_PDCD4)]=0
TPM_PDCD4 <- TPM_PDCD4[, -caret::nearZeroVar(TPM_PDCD4, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PDCD4)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PDCD4), 0.8 * nrow(TPM_PDCD4))
TPM_PDCD4_train <- TPM_PDCD4[train_row,]
TPM_PDCD4_test <- TPM_PDCD4[-train_row,]

dim(TPM_PDCD4_train)
dim(TPM_PDCD4_test)

write.table(TPM_PDCD4_train,"~/Analysis/RBP_validation/split_datasets/TPM_PDCD4_train",sep=";")
write.table(TPM_PDCD4_test,"~/Analysis/RBP_validation/split_datasets/TPM_PDCD4_test",sep=";")

# TPM_PDCD4_train <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_PDCD4_train",sep=";")
# TPM_PDCD4_test <- read.delim("~/Analysis/RBP_validation/split_datasets/TPM_PDCD4_test",sep=";")



## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)

start_time <- Sys.time()
xgb_TPM_PDCD4 <- train(TPM~.,
                        data=TPM_PDCD4_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=60,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PDCD4)

saveRDS(xgb_TPM_PDCD4, "~/Analysis/RBP_validation/models/xgb_TPM_PDCD4_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/RBP_validation/models/xgb_TPM_PDCD4_03_05_2022.RDS")

xgb_TPM_PDCD4_predict_test <- predict(xgb_TPM_PDCD4,TPM_PDCD4_test)
xgb_TPM_PDCD4_predict_test <- data.frame(xgb_TPM_PDCD4_predict_test)
print(lm(xgb_TPM_PDCD4_predict_test$xgb_TPM_PDCD4_predict_test~TPM_PDCD4_test$TPM))
print(cor(xgb_TPM_PDCD4_predict_test$xgb_TPM_PDCD4_predict_test,TPM_PDCD4_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_PDCD4,TPM_PDCD4_param)











