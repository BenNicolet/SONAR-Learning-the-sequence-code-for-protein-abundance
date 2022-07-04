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
###_________________________________________ HEK293T ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data for HEK293T")
# TPM_HEK293T <- read.delim("~/Analysis/data/HEK293T_TPM_log10.csv",sep = ";",dec = ",")
# TPM_HEK293T <- merge(TPM_HEK293T,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_HEK293T_param <- ddply(TPM_HEK293T,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_HEK293T_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_HEK293T$TPM <- 10^(TPM_HEK293T$TPM)
# TPM_HEK293T <- data.frame("ID"=TPM_HEK293T$ensembl_gene_id ,"TPM"=TPM_HEK293T$TPM)
# TPM_HEK293T <- ddply(TPM_HEK293T,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_HEK293T <- merge(TPM_HEK293T,TPM_HEK293T_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_HEK293T$TPM <- log10(TPM_HEK293T$TPM)
# 
# TPM_HEK293T[is.na(TPM_HEK293T)]= 0
# TPM_HEK293T <- subset(TPM_HEK293T, TPM_HEK293T$TPM>=-1)
# print(dim(TPM_HEK293T)) # 13729  2682
# 
# rownames(TPM_HEK293T) <- TPM_HEK293T$ID
# TPM_HEK293T$ID <- NULL
# TPM_HEK293T$ensembl_gene_id <- NULL
# TPM_HEK293T$gene_biotype <- NULL
# TPM_HEK293T$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_HEK293T[is.na(TPM_HEK293T)]=0
# TPM_HEK293T <- TPM_HEK293T[, -caret::nearZeroVar(TPM_HEK293T, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_HEK293T)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_HEK293T), 0.8 * nrow(TPM_HEK293T))
# TPM_HEK293T_train <- TPM_HEK293T[train_row,]
# TPM_HEK293T_test <- TPM_HEK293T[-train_row,]
# 
# dim(TPM_HEK293T_train)
# dim(TPM_HEK293T_test)
# 
# write.table(TPM_HEK293T_train,"~/Analysis/split_dataset/TPM_HEK293T_train",sep=";")
# write.table(TPM_HEK293T_test,"~/Analysis/split_dataset/TPM_HEK293T_test",sep=";")

TPM_HEK293T_train <- read.delim("~/Analysis/split_dataset/TPM_HEK293T_train",sep=";")
TPM_HEK293T_test <- read.delim("~/Analysis/split_dataset/TPM_HEK293T_test",sep=";")



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
xgb_TPM_HEK293T <- train(TPM~.,
                        data=TPM_HEK293T_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=90,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_HEK293T)

saveRDS(xgb_TPM_HEK293T, "~/Analysis/models/xgb_TPM_HEK293T_23_05_2022.RDS")

# xgb_TPM_HEK293T <- readRDS( "~/Analysis/models/xgb_TPM_HEK293T_23_05_2022.RDS")

xgb_TPM_HEK293T_predict_test <- predict(xgb_TPM_HEK293T,TPM_HEK293T_test)
xgb_TPM_HEK293T_predict_test <- data.frame(xgb_TPM_HEK293T_predict_test)
print(lm(xgb_TPM_HEK293T_predict_test$xgb_TPM_HEK293T_predict_test~TPM_HEK293T_test$TPM))
print(cor(xgb_TPM_HEK293T_predict_test$xgb_TPM_HEK293T_predict_test,TPM_HEK293T_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_HEK293T_predict_test$xgb_TPM_HEK293T_predict_test,TPM_HEK293T_test$TPM)

# K562_HEK293T_pred_vs_test_set <- cbind.data.frame(xgb_TPM_HEK293T_predict_test$xgb_TPM_HEK293T_predict_test,TPM_HEK293T_test$TPM)
# 
# write.table(K562_HEK293T_pred_vs_test_set,"/home/nicol01b/Analysis/K562_HEK293T_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_HEK293T,TPM_HEK293T_param)













###__________________________________________________________________________________________###
###_________________________________________ HeLa ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for HeLa")
TPM_HeLa <- read.delim("~/Analysis/data/HeLa_TPM_log10.csv",sep = ";",dec = ",")
TPM_HeLa <- merge(TPM_HeLa,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_HeLa_param <- ddply(TPM_HeLa,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_HeLa_param$TPM <- NULL

print("integrating TPM per gene")
TPM_HeLa$TPM <- 10^(TPM_HeLa$TPM)
TPM_HeLa <- data.frame("ID"=TPM_HeLa$ensembl_gene_id ,"TPM"=TPM_HeLa$TPM)
TPM_HeLa <- ddply(TPM_HeLa,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_HeLa <- merge(TPM_HeLa,TPM_HeLa_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_HeLa$TPM <- log10(TPM_HeLa$TPM)

TPM_HeLa[is.na(TPM_HeLa)]= 0
TPM_HeLa <- subset(TPM_HeLa, TPM_HeLa$TPM>=-1)
print(dim(TPM_HeLa)) # 13729  2682

rownames(TPM_HeLa) <- TPM_HeLa$ID
TPM_HeLa$ID <- NULL
TPM_HeLa$ensembl_gene_id <- NULL
TPM_HeLa$gene_biotype <- NULL
TPM_HeLa$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_HeLa[is.na(TPM_HeLa)]=0
TPM_HeLa <- TPM_HeLa[, -caret::nearZeroVar(TPM_HeLa, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_HeLa)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_HeLa), 0.8 * nrow(TPM_HeLa))
TPM_HeLa_train <- TPM_HeLa[train_row,]
TPM_HeLa_test <- TPM_HeLa[-train_row,]

dim(TPM_HeLa_train)
dim(TPM_HeLa_test)

write.table(TPM_HeLa_train,"~/Analysis/split_dataset/TPM_HeLa_train",sep=";")
write.table(TPM_HeLa_test,"~/Analysis/split_dataset/TPM_HeLa_test",sep=";")
# 
# TPM_HeLa_train <- read.delim("~/Analysis/split_dataset/TPM_HeLa_train",sep=";")
# TPM_HeLa_test <- read.delim("~/Analysis/split_dataset/TPM_HeLa_test",sep=";")



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
xgb_TPM_HeLa <- train(TPM~.,
                         data=TPM_HeLa_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_HeLa)

saveRDS(xgb_TPM_HeLa, "~/Analysis/models/xgb_TPM_HeLa_23_05_2022.RDS")

# xgb_TPM_HeLa <- readRDS( "~/Analysis/models/xgb_TPM_HeLa_23_05_2022.RDS")

xgb_TPM_HeLa_predict_test <- predict(xgb_TPM_HeLa,TPM_HeLa_test)
xgb_TPM_HeLa_predict_test <- data.frame(xgb_TPM_HeLa_predict_test)
print(lm(xgb_TPM_HeLa_predict_test$xgb_TPM_HeLa_predict_test~TPM_HeLa_test$TPM))
print(cor(xgb_TPM_HeLa_predict_test$xgb_TPM_HeLa_predict_test,TPM_HeLa_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_HeLa_predict_test$xgb_TPM_HeLa_predict_test,TPM_HeLa_test$TPM)

# K562_HeLa_pred_vs_test_set <- cbind.data.frame(xgb_TPM_HeLa_predict_test$xgb_TPM_HeLa_predict_test,TPM_HeLa_test$TPM)
# 
# write.table(K562_HeLa_pred_vs_test_set,"/home/nicol01b/Analysis/K562_HeLa_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_HeLa,TPM_HeLa_param)














###__________________________________________________________________________________________###
###_________________________________________ K562 ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for K562")
TPM_K562 <- read.delim("~/Analysis/data/K562_TPM_log10.csv",sep = ";",dec = ",")
TPM_K562 <- merge(TPM_K562,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_K562_param <- ddply(TPM_K562,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_K562_param$TPM <- NULL

print("integrating TPM per gene")
TPM_K562$TPM <- 10^(TPM_K562$TPM)
TPM_K562 <- data.frame("ID"=TPM_K562$ensembl_gene_id ,"TPM"=TPM_K562$TPM)
TPM_K562 <- ddply(TPM_K562,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_K562 <- merge(TPM_K562,TPM_K562_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_K562$TPM <- log10(TPM_K562$TPM)

TPM_K562[is.na(TPM_K562)]= 0
TPM_K562 <- subset(TPM_K562, TPM_K562$TPM>=-1)
print(dim(TPM_K562)) # 13729  2682

rownames(TPM_K562) <- TPM_K562$ID
TPM_K562$ID <- NULL
TPM_K562$ensembl_gene_id <- NULL
TPM_K562$gene_biotype <- NULL
TPM_K562$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_K562[is.na(TPM_K562)]=0
TPM_K562 <- TPM_K562[, -caret::nearZeroVar(TPM_K562, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_K562)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_K562), 0.8 * nrow(TPM_K562))
TPM_K562_train <- TPM_K562[train_row,]
TPM_K562_test <- TPM_K562[-train_row,]

dim(TPM_K562_train)
dim(TPM_K562_test)

write.table(TPM_K562_train,"~/Analysis/split_dataset/TPM_K562_train",sep=";")
write.table(TPM_K562_test,"~/Analysis/split_dataset/TPM_K562_test",sep=";")

# TPM_K562_train <- read.delim("~/Analysis/split_dataset/TPM_K562_train",sep=";")
# TPM_K562_test <- read.delim("~/Analysis/split_dataset/TPM_K562_test",sep=";")



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
xgb_TPM_K562 <- train(TPM~.,
                         data=TPM_K562_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_K562)

saveRDS(xgb_TPM_K562, "~/Analysis/models/xgb_TPM_K562_23_05_2022.RDS")

# xgb_TPM_K562 <- readRDS( "~/Analysis/models/xgb_TPM_K562_23_05_2022.RDS")

xgb_TPM_K562_predict_test <- predict(xgb_TPM_K562,TPM_K562_test)
xgb_TPM_K562_predict_test <- data.frame(xgb_TPM_K562_predict_test)
print(lm(xgb_TPM_K562_predict_test$xgb_TPM_K562_predict_test~TPM_K562_test$TPM))
print(cor(xgb_TPM_K562_predict_test$xgb_TPM_K562_predict_test,TPM_K562_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_K562_predict_test$xgb_TPM_K562_predict_test,TPM_K562_test$TPM)

# K562_K562_pred_vs_test_set <- cbind.data.frame(xgb_TPM_K562_predict_test$xgb_TPM_K562_predict_test,TPM_K562_test$TPM)
# 
# write.table(K562_K562_pred_vs_test_set,"/home/nicol01b/Analysis/K562_K562_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_K562,TPM_K562_param)
























