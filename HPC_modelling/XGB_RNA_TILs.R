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
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "https://apr2018.archive.ensembl.org")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)
#
#
# ## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)




###__________________________________________________________________________________________###
###_________________________________________ EM_A ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for EM_A")
TPM_EM_A <- read.delim("~/Analysis/data/TILs/EM_A.csv",sep = ";",dec = ",")
TPM_EM_A <- merge(TPM_EM_A,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_EM_A_param <- ddply(TPM_EM_A,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_EM_A_param$TPM <- NULL

print("integrating TPM per gene")
TPM_EM_A$TPM <- 10^(TPM_EM_A$TPM)
TPM_EM_A <- data.frame("ID"=TPM_EM_A$ensembl_gene_id ,"TPM"=TPM_EM_A$TPM)
TPM_EM_A <- ddply(TPM_EM_A,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_EM_A <- merge(TPM_EM_A,TPM_EM_A_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_EM_A$TPM <- log10(TPM_EM_A$TPM)

TPM_EM_A[is.na(TPM_EM_A)]= 0
TPM_EM_A <- subset(TPM_EM_A, TPM_EM_A$TPM>=-1)
print(dim(TPM_EM_A)) # 13729  2682

rownames(TPM_EM_A) <- TPM_EM_A$ID
TPM_EM_A$ID <- NULL
TPM_EM_A$ensembl_gene_id <- NULL
TPM_EM_A$gene_biotype <- NULL
TPM_EM_A$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_EM_A[is.na(TPM_EM_A)]=0
TPM_EM_A <- TPM_EM_A[, -caret::nearZeroVar(TPM_EM_A, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_EM_A)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_EM_A), 0.8 * nrow(TPM_EM_A))
TPM_EM_A_train <- TPM_EM_A[train_row,]
TPM_EM_A_test <- TPM_EM_A[-train_row,]

dim(TPM_EM_A_train)
dim(TPM_EM_A_test)

write.table(TPM_EM_A_train,"~/Analysis/split_dataset/TPM_EM_A_train",sep=";")
write.table(TPM_EM_A_test,"~/Analysis/split_dataset/TPM_EM_A_test",sep=";")

# TPM_EM_A_train <- read.delim("~/Analysis/split_dataset/TPM_EM_A_train",sep=";")
# TPM_EM_A_test <- read.delim("~/Analysis/split_dataset/TPM_EM_A_test",sep=";")



## Training the XGB model ##
print("modeling EM_A")
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
xgb_TPM_EM_A <- train(TPM~.,
                         data=TPM_EM_A_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid,
                         na.action = na.omit,
                         nthread=80,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_EM_A)

saveRDS(xgb_TPM_EM_A, "~/Analysis/models/xgb_TPM_EM_A_24_06_2022.RDS")
# xgb_TPM_EM_A <- readRDS( "~/Analysis/models/xgb_TPM_EM_A_24_06_2022.RDS")

xgb_TPM_EM_A_predict_test <- predict(xgb_TPM_EM_A,TPM_EM_A_test)
xgb_TPM_EM_A_predict_test <- data.frame(xgb_TPM_EM_A_predict_test)

print(cor(xgb_TPM_EM_A_predict_test$xgb_TPM_EM_A_predict_test,TPM_EM_A_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_EM_A_predict_test$xgb_TPM_EM_A_predict_test,TPM_EM_A_test$TPM)

# K562_EM_A_pred_vs_test_set <- cbind.data.frame(xgb_TPM_EM_A_predict_test$xgb_TPM_EM_A_predict_test,TPM_EM_A_test$TPM)
# 
# write.table(K562_EM_A_pred_vs_test_set,"/home/nicol01b/Analysis/K562_EM_A_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_EM_A,TPM_EM_A_param)


















###__________________________________________________________________________________________###
###_________________________________________ EM_B ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for EM_B")
TPM_EM_B <- read.delim("~/Analysis/data/TILs/EM_B.csv",sep = ";",dec = ",")
TPM_EM_B <- merge(TPM_EM_B,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_EM_B_param <- ddply(TPM_EM_B,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_EM_B_param$TPM <- NULL

print("integrating TPM per gene")
TPM_EM_B$TPM <- 10^(TPM_EM_B$TPM)
TPM_EM_B <- data.frame("ID"=TPM_EM_B$ensembl_gene_id ,"TPM"=TPM_EM_B$TPM)
TPM_EM_B <- ddply(TPM_EM_B,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_EM_B <- merge(TPM_EM_B,TPM_EM_B_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_EM_B$TPM <- log10(TPM_EM_B$TPM)

TPM_EM_B[is.na(TPM_EM_B)]= 0
TPM_EM_B <- subset(TPM_EM_B, TPM_EM_B$TPM>=-1)
print(dim(TPM_EM_B)) # 13729  2682

rownames(TPM_EM_B) <- TPM_EM_B$ID
TPM_EM_B$ID <- NULL
TPM_EM_B$ensembl_gene_id <- NULL
TPM_EM_B$gene_biotype <- NULL
TPM_EM_B$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_EM_B[is.na(TPM_EM_B)]=0
TPM_EM_B <- TPM_EM_B[, -caret::nearZeroVar(TPM_EM_B, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_EM_B)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_EM_B), 0.8 * nrow(TPM_EM_B))
TPM_EM_B_train <- TPM_EM_B[train_row,]
TPM_EM_B_test <- TPM_EM_B[-train_row,]

dim(TPM_EM_B_train)
dim(TPM_EM_B_test)

write.table(TPM_EM_B_train,"~/Analysis/split_dataset/TPM_EM_B_train",sep=";")
write.table(TPM_EM_B_test,"~/Analysis/split_dataset/TPM_EM_B_test",sep=";")

# TPM_EM_B_train <- read.delim("~/Analysis/split_dataset/TPM_EM_B_train",sep=";")
# TPM_EM_B_test <- read.delim("~/Analysis/split_dataset/TPM_EM_B_test",sep=";")



## Training the XGB model ##
print("modeling EM_B")
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
xgb_TPM_EM_B <- train(TPM~.,
                      data=TPM_EM_B_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_EM_B)

saveRDS(xgb_TPM_EM_B, "~/Analysis/models/xgb_TPM_EM_B_24_06_2022.RDS")
# xgb_TPM_EM_B <- readRDS( "~/Analysis/models/xgb_TPM_EM_B_24_06_2022.RDS")

xgb_TPM_EM_B_predict_test <- predict(xgb_TPM_EM_B,TPM_EM_B_test)
xgb_TPM_EM_B_predict_test <- data.frame(xgb_TPM_EM_B_predict_test)

print(cor(xgb_TPM_EM_B_predict_test$xgb_TPM_EM_B_predict_test,TPM_EM_B_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_EM_B_predict_test$xgb_TPM_EM_B_predict_test,TPM_EM_B_test$TPM)

# K562_EM_B_pred_vs_test_set <- cbind.data.frame(xgb_TPM_EM_B_predict_test$xgb_TPM_EM_B_predict_test,TPM_EM_B_test$TPM)
# 
# write.table(K562_EM_B_pred_vs_test_set,"/home/nicol01b/Analysis/K562_EM_B_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_EM_B,TPM_EM_B_param)

















###__________________________________________________________________________________________###
###_________________________________________ EM_C ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for EM_C")
TPM_EM_C <- read.delim("~/Analysis/data/TILs/EM_C.csv",sep = ";",dec = ",")
TPM_EM_C <- merge(TPM_EM_C,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_EM_C_param <- ddply(TPM_EM_C,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_EM_C_param$TPM <- NULL

print("integrating TPM per gene")
TPM_EM_C$TPM <- 10^(TPM_EM_C$TPM)
TPM_EM_C <- data.frame("ID"=TPM_EM_C$ensembl_gene_id ,"TPM"=TPM_EM_C$TPM)
TPM_EM_C <- ddply(TPM_EM_C,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_EM_C <- merge(TPM_EM_C,TPM_EM_C_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_EM_C$TPM <- log10(TPM_EM_C$TPM)

TPM_EM_C[is.na(TPM_EM_C)]= 0
TPM_EM_C <- subset(TPM_EM_C, TPM_EM_C$TPM>=-1)
print(dim(TPM_EM_C)) # 13729  2682

rownames(TPM_EM_C) <- TPM_EM_C$ID
TPM_EM_C$ID <- NULL
TPM_EM_C$ensembl_gene_id <- NULL
TPM_EM_C$gene_biotype <- NULL
TPM_EM_C$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_EM_C[is.na(TPM_EM_C)]=0
TPM_EM_C <- TPM_EM_C[, -caret::nearZeroVar(TPM_EM_C, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_EM_C)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_EM_C), 0.8 * nrow(TPM_EM_C))
TPM_EM_C_train <- TPM_EM_C[train_row,]
TPM_EM_C_test <- TPM_EM_C[-train_row,]

dim(TPM_EM_C_train)
dim(TPM_EM_C_test)

write.table(TPM_EM_C_train,"~/Analysis/split_dataset/TPM_EM_C_train",sep=";")
write.table(TPM_EM_C_test,"~/Analysis/split_dataset/TPM_EM_C_test",sep=";")

# TPM_EM_C_train <- read.delim("~/Analysis/split_dataset/TPM_EM_C_train",sep=";")
# TPM_EM_C_test <- read.delim("~/Analysis/split_dataset/TPM_EM_C_test",sep=";")



## Training the XGB model ##
print("modeling EM_C")
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
xgb_TPM_EM_C <- train(TPM~.,
                      data=TPM_EM_C_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_EM_C)

saveRDS(xgb_TPM_EM_C, "~/Analysis/models/xgb_TPM_EM_C_24_06_2022.RDS")
# xgb_TPM_EM_C <- readRDS( "~/Analysis/models/xgb_TPM_EM_C_24_06_2022.RDS")

xgb_TPM_EM_C_predict_test <- predict(xgb_TPM_EM_C,TPM_EM_C_test)
xgb_TPM_EM_C_predict_test <- data.frame(xgb_TPM_EM_C_predict_test)

print(cor(xgb_TPM_EM_C_predict_test$xgb_TPM_EM_C_predict_test,TPM_EM_C_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_EM_C_predict_test$xgb_TPM_EM_C_predict_test,TPM_EM_C_test$TPM)

# K562_EM_C_pred_vs_test_set <- cbind.data.frame(xgb_TPM_EM_C_predict_test$xgb_TPM_EM_C_predict_test,TPM_EM_C_test$TPM)
# 
# write.table(K562_EM_C_pred_vs_test_set,"/home/nicol01b/Analysis/K562_EM_C_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_EM_C,TPM_EM_C_param)














###__________________________________________________________________________________________###
###_________________________________________ EM_D ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for EM_D")
TPM_EM_D <- read.delim("~/Analysis/data/TILs/EM_D.csv",sep = ";",dec = ",")
TPM_EM_D <- merge(TPM_EM_D,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_EM_D_param <- ddply(TPM_EM_D,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_EM_D_param$TPM <- NULL

print("integrating TPM per gene")
TPM_EM_D$TPM <- 10^(TPM_EM_D$TPM)
TPM_EM_D <- data.frame("ID"=TPM_EM_D$ensembl_gene_id ,"TPM"=TPM_EM_D$TPM)
TPM_EM_D <- ddply(TPM_EM_D,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_EM_D <- merge(TPM_EM_D,TPM_EM_D_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_EM_D$TPM <- log10(TPM_EM_D$TPM)

TPM_EM_D[is.na(TPM_EM_D)]= 0
TPM_EM_D <- subset(TPM_EM_D, TPM_EM_D$TPM>=-1)
print(dim(TPM_EM_D)) # 13729  2682

rownames(TPM_EM_D) <- TPM_EM_D$ID
TPM_EM_D$ID <- NULL
TPM_EM_D$ensembl_gene_id <- NULL
TPM_EM_D$gene_biotype <- NULL
TPM_EM_D$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_EM_D[is.na(TPM_EM_D)]=0
TPM_EM_D <- TPM_EM_D[, -caret::nearZeroVar(TPM_EM_D, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_EM_D)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_EM_D), 0.8 * nrow(TPM_EM_D))
TPM_EM_D_train <- TPM_EM_D[train_row,]
TPM_EM_D_test <- TPM_EM_D[-train_row,]

dim(TPM_EM_D_train)
dim(TPM_EM_D_test)

write.table(TPM_EM_D_train,"~/Analysis/split_dataset/TPM_EM_D_train",sep=";")
write.table(TPM_EM_D_test,"~/Analysis/split_dataset/TPM_EM_D_test",sep=";")

# TPM_EM_D_train <- read.delim("~/Analysis/split_dataset/TPM_EM_D_train",sep=";")
# TPM_EM_D_test <- read.delim("~/Analysis/split_dataset/TPM_EM_D_test",sep=";")



## Training the XGB model ##
print("modeling EM_D")
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
xgb_TPM_EM_D <- train(TPM~.,
                      data=TPM_EM_D_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_EM_D)

saveRDS(xgb_TPM_EM_D, "~/Analysis/models/xgb_TPM_EM_D_24_06_2022.RDS")
# xgb_TPM_EM_D <- readRDS( "~/Analysis/models/xgb_TPM_EM_D_24_06_2022.RDS")

xgb_TPM_EM_D_predict_test <- predict(xgb_TPM_EM_D,TPM_EM_D_test)
xgb_TPM_EM_D_predict_test <- data.frame(xgb_TPM_EM_D_predict_test)

print(cor(xgb_TPM_EM_D_predict_test$xgb_TPM_EM_D_predict_test,TPM_EM_D_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_EM_D_predict_test$xgb_TPM_EM_D_predict_test,TPM_EM_D_test$TPM)

# K562_EM_D_pred_vs_test_set <- cbind.data.frame(xgb_TPM_EM_D_predict_test$xgb_TPM_EM_D_predict_test,TPM_EM_D_test$TPM)
# 
# write.table(K562_EM_D_pred_vs_test_set,"/home/nicol01b/Analysis/K562_EM_D_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_EM_D,TPM_EM_D_param)


















##___________________________________________________________________________________________________________________________##













###__________________________________________________________________________________________###
###_________________________________________ PD1_hi_A ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_hi_A")
TPM_PD1_hi_A <- read.delim("~/Analysis/data/TILs/PD1_hi_A.csv",sep = ";",dec = ",")
TPM_PD1_hi_A <- merge(TPM_PD1_hi_A,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_hi_A_param <- ddply(TPM_PD1_hi_A,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_hi_A_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_hi_A$TPM <- 10^(TPM_PD1_hi_A$TPM)
TPM_PD1_hi_A <- data.frame("ID"=TPM_PD1_hi_A$ensembl_gene_id ,"TPM"=TPM_PD1_hi_A$TPM)
TPM_PD1_hi_A <- ddply(TPM_PD1_hi_A,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_hi_A <- merge(TPM_PD1_hi_A,TPM_PD1_hi_A_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_hi_A$TPM <- log10(TPM_PD1_hi_A$TPM)

TPM_PD1_hi_A[is.na(TPM_PD1_hi_A)]= 0
TPM_PD1_hi_A <- subset(TPM_PD1_hi_A, TPM_PD1_hi_A$TPM>=-1)
print(dim(TPM_PD1_hi_A)) # 13729  2682

rownames(TPM_PD1_hi_A) <- TPM_PD1_hi_A$ID
TPM_PD1_hi_A$ID <- NULL
TPM_PD1_hi_A$ensembl_gene_id <- NULL
TPM_PD1_hi_A$gene_biotype <- NULL
TPM_PD1_hi_A$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_hi_A[is.na(TPM_PD1_hi_A)]=0
TPM_PD1_hi_A <- TPM_PD1_hi_A[, -caret::nearZeroVar(TPM_PD1_hi_A, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_hi_A)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_hi_A), 0.8 * nrow(TPM_PD1_hi_A))
TPM_PD1_hi_A_train <- TPM_PD1_hi_A[train_row,]
TPM_PD1_hi_A_test <- TPM_PD1_hi_A[-train_row,]

dim(TPM_PD1_hi_A_train)
dim(TPM_PD1_hi_A_test)

write.table(TPM_PD1_hi_A_train,"~/Analysis/split_dataset/TPM_PD1_hi_A_train",sep=";")
write.table(TPM_PD1_hi_A_test,"~/Analysis/split_dataset/TPM_PD1_hi_A_test",sep=";")

# TPM_PD1_hi_A_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_A_train",sep=";")
# TPM_PD1_hi_A_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_A_test",sep=";")



## Training the XGB model ##
print("modeling PD1_hi_A")
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
xgb_TPM_PD1_hi_A <- train(TPM~.,
                      data=TPM_PD1_hi_A_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_hi_A)

saveRDS(xgb_TPM_PD1_hi_A, "~/Analysis/models/xgb_TPM_PD1_hi_A_24_06_2022.RDS")
# xgb_TPM_PD1_hi_A <- readRDS( "~/Analysis/models/xgb_TPM_PD1_hi_A_24_06_2022.RDS")

xgb_TPM_PD1_hi_A_predict_test <- predict(xgb_TPM_PD1_hi_A,TPM_PD1_hi_A_test)
xgb_TPM_PD1_hi_A_predict_test <- data.frame(xgb_TPM_PD1_hi_A_predict_test)

print(cor(xgb_TPM_PD1_hi_A_predict_test$xgb_TPM_PD1_hi_A_predict_test,TPM_PD1_hi_A_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_hi_A_predict_test$xgb_TPM_PD1_hi_A_predict_test,TPM_PD1_hi_A_test$TPM)

# K562_PD1_hi_A_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_hi_A_predict_test$xgb_TPM_PD1_hi_A_predict_test,TPM_PD1_hi_A_test$TPM)
# 
# write.table(K562_PD1_hi_A_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_hi_A_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_hi_A,TPM_PD1_hi_A_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_hi_B ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_hi_B")
TPM_PD1_hi_B <- read.delim("~/Analysis/data/TILs/PD1_hi_B.csv",sep = ";",dec = ",")
TPM_PD1_hi_B <- merge(TPM_PD1_hi_B,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_hi_B_param <- ddply(TPM_PD1_hi_B,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_hi_B_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_hi_B$TPM <- 10^(TPM_PD1_hi_B$TPM)
TPM_PD1_hi_B <- data.frame("ID"=TPM_PD1_hi_B$ensembl_gene_id ,"TPM"=TPM_PD1_hi_B$TPM)
TPM_PD1_hi_B <- ddply(TPM_PD1_hi_B,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_hi_B <- merge(TPM_PD1_hi_B,TPM_PD1_hi_B_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_hi_B$TPM <- log10(TPM_PD1_hi_B$TPM)

TPM_PD1_hi_B[is.na(TPM_PD1_hi_B)]= 0
TPM_PD1_hi_B <- subset(TPM_PD1_hi_B, TPM_PD1_hi_B$TPM>=-1)
print(dim(TPM_PD1_hi_B)) # 13729  2682

rownames(TPM_PD1_hi_B) <- TPM_PD1_hi_B$ID
TPM_PD1_hi_B$ID <- NULL
TPM_PD1_hi_B$ensembl_gene_id <- NULL
TPM_PD1_hi_B$gene_biotype <- NULL
TPM_PD1_hi_B$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_hi_B[is.na(TPM_PD1_hi_B)]=0
TPM_PD1_hi_B <- TPM_PD1_hi_B[, -caret::nearZeroVar(TPM_PD1_hi_B, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_hi_B)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_hi_B), 0.8 * nrow(TPM_PD1_hi_B))
TPM_PD1_hi_B_train <- TPM_PD1_hi_B[train_row,]
TPM_PD1_hi_B_test <- TPM_PD1_hi_B[-train_row,]

dim(TPM_PD1_hi_B_train)
dim(TPM_PD1_hi_B_test)

write.table(TPM_PD1_hi_B_train,"~/Analysis/split_dataset/TPM_PD1_hi_B_train",sep=";")
write.table(TPM_PD1_hi_B_test,"~/Analysis/split_dataset/TPM_PD1_hi_B_test",sep=";")

# TPM_PD1_hi_B_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_B_train",sep=";")
# TPM_PD1_hi_B_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_B_test",sep=";")



## Training the XGB model ##
print("modeling PD1_hi_B")
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
xgb_TPM_PD1_hi_B <- train(TPM~.,
                      data=TPM_PD1_hi_B_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_hi_B)

saveRDS(xgb_TPM_PD1_hi_B, "~/Analysis/models/xgb_TPM_PD1_hi_B_24_06_2022.RDS")
# xgb_TPM_PD1_hi_B <- readRDS( "~/Analysis/models/xgb_TPM_PD1_hi_B_24_06_2022.RDS")

xgb_TPM_PD1_hi_B_predict_test <- predict(xgb_TPM_PD1_hi_B,TPM_PD1_hi_B_test)
xgb_TPM_PD1_hi_B_predict_test <- data.frame(xgb_TPM_PD1_hi_B_predict_test)

print(cor(xgb_TPM_PD1_hi_B_predict_test$xgb_TPM_PD1_hi_B_predict_test,TPM_PD1_hi_B_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_hi_B_predict_test$xgb_TPM_PD1_hi_B_predict_test,TPM_PD1_hi_B_test$TPM)

# K562_PD1_hi_B_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_hi_B_predict_test$xgb_TPM_PD1_hi_B_predict_test,TPM_PD1_hi_B_test$TPM)
# 
# write.table(K562_PD1_hi_B_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_hi_B_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_hi_B,TPM_PD1_hi_B_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_hi_C ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_hi_C")
TPM_PD1_hi_C <- read.delim("~/Analysis/data/TILs/PD1_hi_C.csv",sep = ";",dec = ",")
TPM_PD1_hi_C <- merge(TPM_PD1_hi_C,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_hi_C_param <- ddply(TPM_PD1_hi_C,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_hi_C_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_hi_C$TPM <- 10^(TPM_PD1_hi_C$TPM)
TPM_PD1_hi_C <- data.frame("ID"=TPM_PD1_hi_C$ensembl_gene_id ,"TPM"=TPM_PD1_hi_C$TPM)
TPM_PD1_hi_C <- ddply(TPM_PD1_hi_C,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_hi_C <- merge(TPM_PD1_hi_C,TPM_PD1_hi_C_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_hi_C$TPM <- log10(TPM_PD1_hi_C$TPM)

TPM_PD1_hi_C[is.na(TPM_PD1_hi_C)]= 0
TPM_PD1_hi_C <- subset(TPM_PD1_hi_C, TPM_PD1_hi_C$TPM>=-1)
print(dim(TPM_PD1_hi_C)) # 13729  2682

rownames(TPM_PD1_hi_C) <- TPM_PD1_hi_C$ID
TPM_PD1_hi_C$ID <- NULL
TPM_PD1_hi_C$ensembl_gene_id <- NULL
TPM_PD1_hi_C$gene_biotype <- NULL
TPM_PD1_hi_C$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_hi_C[is.na(TPM_PD1_hi_C)]=0
TPM_PD1_hi_C <- TPM_PD1_hi_C[, -caret::nearZeroVar(TPM_PD1_hi_C, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_hi_C)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_hi_C), 0.8 * nrow(TPM_PD1_hi_C))
TPM_PD1_hi_C_train <- TPM_PD1_hi_C[train_row,]
TPM_PD1_hi_C_test <- TPM_PD1_hi_C[-train_row,]

dim(TPM_PD1_hi_C_train)
dim(TPM_PD1_hi_C_test)

write.table(TPM_PD1_hi_C_train,"~/Analysis/split_dataset/TPM_PD1_hi_C_train",sep=";")
write.table(TPM_PD1_hi_C_test,"~/Analysis/split_dataset/TPM_PD1_hi_C_test",sep=";")

# TPM_PD1_hi_C_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_C_train",sep=";")
# TPM_PD1_hi_C_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_C_test",sep=";")



## Training the XGB model ##
print("modeling PD1_hi_C")
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
xgb_TPM_PD1_hi_C <- train(TPM~.,
                      data=TPM_PD1_hi_C_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_hi_C)

saveRDS(xgb_TPM_PD1_hi_C, "~/Analysis/models/xgb_TPM_PD1_hi_C_24_06_2022.RDS")
# xgb_TPM_PD1_hi_C <- readRDS( "~/Analysis/models/xgb_TPM_PD1_hi_C_24_06_2022.RDS")

xgb_TPM_PD1_hi_C_predict_test <- predict(xgb_TPM_PD1_hi_C,TPM_PD1_hi_C_test)
xgb_TPM_PD1_hi_C_predict_test <- data.frame(xgb_TPM_PD1_hi_C_predict_test)

print(cor(xgb_TPM_PD1_hi_C_predict_test$xgb_TPM_PD1_hi_C_predict_test,TPM_PD1_hi_C_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_hi_C_predict_test$xgb_TPM_PD1_hi_C_predict_test,TPM_PD1_hi_C_test$TPM)

# K562_PD1_hi_C_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_hi_C_predict_test$xgb_TPM_PD1_hi_C_predict_test,TPM_PD1_hi_C_test$TPM)
# 
# write.table(K562_PD1_hi_C_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_hi_C_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_hi_C,TPM_PD1_hi_C_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_hi_D ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_hi_D")
TPM_PD1_hi_D <- read.delim("~/Analysis/data/TILs/PD1_hi_D.csv",sep = ";",dec = ",")
TPM_PD1_hi_D <- merge(TPM_PD1_hi_D,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_hi_D_param <- ddply(TPM_PD1_hi_D,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_hi_D_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_hi_D$TPM <- 10^(TPM_PD1_hi_D$TPM)
TPM_PD1_hi_D <- data.frame("ID"=TPM_PD1_hi_D$ensembl_gene_id ,"TPM"=TPM_PD1_hi_D$TPM)
TPM_PD1_hi_D <- ddply(TPM_PD1_hi_D,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_hi_D <- merge(TPM_PD1_hi_D,TPM_PD1_hi_D_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_hi_D$TPM <- log10(TPM_PD1_hi_D$TPM)

TPM_PD1_hi_D[is.na(TPM_PD1_hi_D)]= 0
TPM_PD1_hi_D <- subset(TPM_PD1_hi_D, TPM_PD1_hi_D$TPM>=-1)
print(dim(TPM_PD1_hi_D)) # 13729  2682

rownames(TPM_PD1_hi_D) <- TPM_PD1_hi_D$ID
TPM_PD1_hi_D$ID <- NULL
TPM_PD1_hi_D$ensembl_gene_id <- NULL
TPM_PD1_hi_D$gene_biotype <- NULL
TPM_PD1_hi_D$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_hi_D[is.na(TPM_PD1_hi_D)]=0
TPM_PD1_hi_D <- TPM_PD1_hi_D[, -caret::nearZeroVar(TPM_PD1_hi_D, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_hi_D)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_hi_D), 0.8 * nrow(TPM_PD1_hi_D))
TPM_PD1_hi_D_train <- TPM_PD1_hi_D[train_row,]
TPM_PD1_hi_D_test <- TPM_PD1_hi_D[-train_row,]

dim(TPM_PD1_hi_D_train)
dim(TPM_PD1_hi_D_test)

write.table(TPM_PD1_hi_D_train,"~/Analysis/split_dataset/TPM_PD1_hi_D_train",sep=";")
write.table(TPM_PD1_hi_D_test,"~/Analysis/split_dataset/TPM_PD1_hi_D_test",sep=";")

# TPM_PD1_hi_D_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_D_train",sep=";")
# TPM_PD1_hi_D_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_hi_D_test",sep=";")



## Training the XGB model ##
print("modeling PD1_hi_D")
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
xgb_TPM_PD1_hi_D <- train(TPM~.,
                      data=TPM_PD1_hi_D_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_hi_D)

saveRDS(xgb_TPM_PD1_hi_D, "~/Analysis/models/xgb_TPM_PD1_hi_D_24_06_2022.RDS")
# xgb_TPM_PD1_hi_D <- readRDS( "~/Analysis/models/xgb_TPM_PD1_hi_D_24_06_2022.RDS")

xgb_TPM_PD1_hi_D_predict_test <- predict(xgb_TPM_PD1_hi_D,TPM_PD1_hi_D_test)
xgb_TPM_PD1_hi_D_predict_test <- data.frame(xgb_TPM_PD1_hi_D_predict_test)

print(cor(xgb_TPM_PD1_hi_D_predict_test$xgb_TPM_PD1_hi_D_predict_test,TPM_PD1_hi_D_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_hi_D_predict_test$xgb_TPM_PD1_hi_D_predict_test,TPM_PD1_hi_D_test$TPM)

# K562_PD1_hi_D_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_hi_D_predict_test$xgb_TPM_PD1_hi_D_predict_test,TPM_PD1_hi_D_test$TPM)
# 
# write.table(K562_PD1_hi_D_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_hi_D_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_hi_D,TPM_PD1_hi_D_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_int_A ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_int_A")
TPM_PD1_int_A <- read.delim("~/Analysis/data/TILs/PD1_int_A.csv",sep = ";",dec = ",")
TPM_PD1_int_A <- merge(TPM_PD1_int_A,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_int_A_param <- ddply(TPM_PD1_int_A,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_int_A_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_int_A$TPM <- 10^(TPM_PD1_int_A$TPM)
TPM_PD1_int_A <- data.frame("ID"=TPM_PD1_int_A$ensembl_gene_id ,"TPM"=TPM_PD1_int_A$TPM)
TPM_PD1_int_A <- ddply(TPM_PD1_int_A,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_int_A <- merge(TPM_PD1_int_A,TPM_PD1_int_A_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_int_A$TPM <- log10(TPM_PD1_int_A$TPM)

TPM_PD1_int_A[is.na(TPM_PD1_int_A)]= 0
TPM_PD1_int_A <- subset(TPM_PD1_int_A, TPM_PD1_int_A$TPM>=-1)
print(dim(TPM_PD1_int_A)) # 13729  2682

rownames(TPM_PD1_int_A) <- TPM_PD1_int_A$ID
TPM_PD1_int_A$ID <- NULL
TPM_PD1_int_A$ensembl_gene_id <- NULL
TPM_PD1_int_A$gene_biotype <- NULL
TPM_PD1_int_A$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_int_A[is.na(TPM_PD1_int_A)]=0
TPM_PD1_int_A <- TPM_PD1_int_A[, -caret::nearZeroVar(TPM_PD1_int_A, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_int_A)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_int_A), 0.8 * nrow(TPM_PD1_int_A))
TPM_PD1_int_A_train <- TPM_PD1_int_A[train_row,]
TPM_PD1_int_A_test <- TPM_PD1_int_A[-train_row,]

dim(TPM_PD1_int_A_train)
dim(TPM_PD1_int_A_test)

write.table(TPM_PD1_int_A_train,"~/Analysis/split_dataset/TPM_PD1_int_A_train",sep=";")
write.table(TPM_PD1_int_A_test,"~/Analysis/split_dataset/TPM_PD1_int_A_test",sep=";")

# TPM_PD1_int_A_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_A_train",sep=";")
# TPM_PD1_int_A_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_A_test",sep=";")



## Training the XGB model ##
print("modeling PD1_int_A")
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
xgb_TPM_PD1_int_A <- train(TPM~.,
                      data=TPM_PD1_int_A_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_int_A)

saveRDS(xgb_TPM_PD1_int_A, "~/Analysis/models/xgb_TPM_PD1_int_A_24_06_2022.RDS")
# xgb_TPM_PD1_int_A <- readRDS( "~/Analysis/models/xgb_TPM_PD1_int_A_24_06_2022.RDS")

xgb_TPM_PD1_int_A_predict_test <- predict(xgb_TPM_PD1_int_A,TPM_PD1_int_A_test)
xgb_TPM_PD1_int_A_predict_test <- data.frame(xgb_TPM_PD1_int_A_predict_test)

print(cor(xgb_TPM_PD1_int_A_predict_test$xgb_TPM_PD1_int_A_predict_test,TPM_PD1_int_A_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_int_A_predict_test$xgb_TPM_PD1_int_A_predict_test,TPM_PD1_int_A_test$TPM)

# K562_PD1_int_A_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_int_A_predict_test$xgb_TPM_PD1_int_A_predict_test,TPM_PD1_int_A_test$TPM)
# 
# write.table(K562_PD1_int_A_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_int_A_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_int_A,TPM_PD1_int_A_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_int_B ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_int_B")
TPM_PD1_int_B <- read.delim("~/Analysis/data/TILs/PD1_int_B.csv",sep = ";",dec = ",")
TPM_PD1_int_B <- merge(TPM_PD1_int_B,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_int_B_param <- ddply(TPM_PD1_int_B,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_int_B_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_int_B$TPM <- 10^(TPM_PD1_int_B$TPM)
TPM_PD1_int_B <- data.frame("ID"=TPM_PD1_int_B$ensembl_gene_id ,"TPM"=TPM_PD1_int_B$TPM)
TPM_PD1_int_B <- ddply(TPM_PD1_int_B,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_int_B <- merge(TPM_PD1_int_B,TPM_PD1_int_B_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_int_B$TPM <- log10(TPM_PD1_int_B$TPM)

TPM_PD1_int_B[is.na(TPM_PD1_int_B)]= 0
TPM_PD1_int_B <- subset(TPM_PD1_int_B, TPM_PD1_int_B$TPM>=-1)
print(dim(TPM_PD1_int_B)) # 13729  2682

rownames(TPM_PD1_int_B) <- TPM_PD1_int_B$ID
TPM_PD1_int_B$ID <- NULL
TPM_PD1_int_B$ensembl_gene_id <- NULL
TPM_PD1_int_B$gene_biotype <- NULL
TPM_PD1_int_B$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_int_B[is.na(TPM_PD1_int_B)]=0
TPM_PD1_int_B <- TPM_PD1_int_B[, -caret::nearZeroVar(TPM_PD1_int_B, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_int_B)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_int_B), 0.8 * nrow(TPM_PD1_int_B))
TPM_PD1_int_B_train <- TPM_PD1_int_B[train_row,]
TPM_PD1_int_B_test <- TPM_PD1_int_B[-train_row,]

dim(TPM_PD1_int_B_train)
dim(TPM_PD1_int_B_test)

write.table(TPM_PD1_int_B_train,"~/Analysis/split_dataset/TPM_PD1_int_B_train",sep=";")
write.table(TPM_PD1_int_B_test,"~/Analysis/split_dataset/TPM_PD1_int_B_test",sep=";")

# TPM_PD1_int_B_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_B_train",sep=";")
# TPM_PD1_int_B_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_B_test",sep=";")



## Training the XGB model ##
print("modeling PD1_int_B")
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
xgb_TPM_PD1_int_B <- train(TPM~.,
                      data=TPM_PD1_int_B_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_int_B)

saveRDS(xgb_TPM_PD1_int_B, "~/Analysis/models/xgb_TPM_PD1_int_B_24_06_2022.RDS")
# xgb_TPM_PD1_int_B <- readRDS( "~/Analysis/models/xgb_TPM_PD1_int_B_24_06_2022.RDS")

xgb_TPM_PD1_int_B_predict_test <- predict(xgb_TPM_PD1_int_B,TPM_PD1_int_B_test)
xgb_TPM_PD1_int_B_predict_test <- data.frame(xgb_TPM_PD1_int_B_predict_test)

print(cor(xgb_TPM_PD1_int_B_predict_test$xgb_TPM_PD1_int_B_predict_test,TPM_PD1_int_B_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_int_B_predict_test$xgb_TPM_PD1_int_B_predict_test,TPM_PD1_int_B_test$TPM)

# K562_PD1_int_B_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_int_B_predict_test$xgb_TPM_PD1_int_B_predict_test,TPM_PD1_int_B_test$TPM)
# 
# write.table(K562_PD1_int_B_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_int_B_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_int_B,TPM_PD1_int_B_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_int_C ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_int_C")
TPM_PD1_int_C <- read.delim("~/Analysis/data/TILs/PD1_int_C.csv",sep = ";",dec = ",")
TPM_PD1_int_C <- merge(TPM_PD1_int_C,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_int_C_param <- ddply(TPM_PD1_int_C,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_int_C_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_int_C$TPM <- 10^(TPM_PD1_int_C$TPM)
TPM_PD1_int_C <- data.frame("ID"=TPM_PD1_int_C$ensembl_gene_id ,"TPM"=TPM_PD1_int_C$TPM)
TPM_PD1_int_C <- ddply(TPM_PD1_int_C,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_int_C <- merge(TPM_PD1_int_C,TPM_PD1_int_C_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_int_C$TPM <- log10(TPM_PD1_int_C$TPM)

TPM_PD1_int_C[is.na(TPM_PD1_int_C)]= 0
TPM_PD1_int_C <- subset(TPM_PD1_int_C, TPM_PD1_int_C$TPM>=-1)
print(dim(TPM_PD1_int_C)) # 13729  2682

rownames(TPM_PD1_int_C) <- TPM_PD1_int_C$ID
TPM_PD1_int_C$ID <- NULL
TPM_PD1_int_C$ensembl_gene_id <- NULL
TPM_PD1_int_C$gene_biotype <- NULL
TPM_PD1_int_C$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_int_C[is.na(TPM_PD1_int_C)]=0
TPM_PD1_int_C <- TPM_PD1_int_C[, -caret::nearZeroVar(TPM_PD1_int_C, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_int_C)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_int_C), 0.8 * nrow(TPM_PD1_int_C))
TPM_PD1_int_C_train <- TPM_PD1_int_C[train_row,]
TPM_PD1_int_C_test <- TPM_PD1_int_C[-train_row,]

dim(TPM_PD1_int_C_train)
dim(TPM_PD1_int_C_test)

write.table(TPM_PD1_int_C_train,"~/Analysis/split_dataset/TPM_PD1_int_C_train",sep=";")
write.table(TPM_PD1_int_C_test,"~/Analysis/split_dataset/TPM_PD1_int_C_test",sep=";")

# TPM_PD1_int_C_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_C_train",sep=";")
# TPM_PD1_int_C_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_C_test",sep=";")



## Training the XGB model ##
print("modeling PD1_int_C")
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
xgb_TPM_PD1_int_C <- train(TPM~.,
                      data=TPM_PD1_int_C_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_int_C)

saveRDS(xgb_TPM_PD1_int_C, "~/Analysis/models/xgb_TPM_PD1_int_C_24_06_2022.RDS")
# xgb_TPM_PD1_int_C <- readRDS( "~/Analysis/models/xgb_TPM_PD1_int_C_24_06_2022.RDS")

xgb_TPM_PD1_int_C_predict_test <- predict(xgb_TPM_PD1_int_C,TPM_PD1_int_C_test)
xgb_TPM_PD1_int_C_predict_test <- data.frame(xgb_TPM_PD1_int_C_predict_test)

print(cor(xgb_TPM_PD1_int_C_predict_test$xgb_TPM_PD1_int_C_predict_test,TPM_PD1_int_C_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_int_C_predict_test$xgb_TPM_PD1_int_C_predict_test,TPM_PD1_int_C_test$TPM)

# K562_PD1_int_C_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_int_C_predict_test$xgb_TPM_PD1_int_C_predict_test,TPM_PD1_int_C_test$TPM)
# 
# write.table(K562_PD1_int_C_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_int_C_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_int_C,TPM_PD1_int_C_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_int_D ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_int_D")
TPM_PD1_int_D <- read.delim("~/Analysis/data/TILs/PD1_int_D.csv",sep = ";",dec = ",")
TPM_PD1_int_D <- merge(TPM_PD1_int_D,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_int_D_param <- ddply(TPM_PD1_int_D,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_int_D_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_int_D$TPM <- 10^(TPM_PD1_int_D$TPM)
TPM_PD1_int_D <- data.frame("ID"=TPM_PD1_int_D$ensembl_gene_id ,"TPM"=TPM_PD1_int_D$TPM)
TPM_PD1_int_D <- ddply(TPM_PD1_int_D,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_int_D <- merge(TPM_PD1_int_D,TPM_PD1_int_D_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_int_D$TPM <- log10(TPM_PD1_int_D$TPM)

TPM_PD1_int_D[is.na(TPM_PD1_int_D)]= 0
TPM_PD1_int_D <- subset(TPM_PD1_int_D, TPM_PD1_int_D$TPM>=-1)
print(dim(TPM_PD1_int_D)) # 13729  2682

rownames(TPM_PD1_int_D) <- TPM_PD1_int_D$ID
TPM_PD1_int_D$ID <- NULL
TPM_PD1_int_D$ensembl_gene_id <- NULL
TPM_PD1_int_D$gene_biotype <- NULL
TPM_PD1_int_D$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_int_D[is.na(TPM_PD1_int_D)]=0
TPM_PD1_int_D <- TPM_PD1_int_D[, -caret::nearZeroVar(TPM_PD1_int_D, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_int_D)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_int_D), 0.8 * nrow(TPM_PD1_int_D))
TPM_PD1_int_D_train <- TPM_PD1_int_D[train_row,]
TPM_PD1_int_D_test <- TPM_PD1_int_D[-train_row,]

dim(TPM_PD1_int_D_train)
dim(TPM_PD1_int_D_test)

write.table(TPM_PD1_int_D_train,"~/Analysis/split_dataset/TPM_PD1_int_D_train",sep=";")
write.table(TPM_PD1_int_D_test,"~/Analysis/split_dataset/TPM_PD1_int_D_test",sep=";")

# TPM_PD1_int_D_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_D_train",sep=";")
# TPM_PD1_int_D_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_int_D_test",sep=";")



## Training the XGB model ##
print("modeling PD1_int_D")
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
xgb_TPM_PD1_int_D <- train(TPM~.,
                      data=TPM_PD1_int_D_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_int_D)

saveRDS(xgb_TPM_PD1_int_D, "~/Analysis/models/xgb_TPM_PD1_int_D_24_06_2022.RDS")
# xgb_TPM_PD1_int_D <- readRDS( "~/Analysis/models/xgb_TPM_PD1_int_D_24_06_2022.RDS")

xgb_TPM_PD1_int_D_predict_test <- predict(xgb_TPM_PD1_int_D,TPM_PD1_int_D_test)
xgb_TPM_PD1_int_D_predict_test <- data.frame(xgb_TPM_PD1_int_D_predict_test)

print(cor(xgb_TPM_PD1_int_D_predict_test$xgb_TPM_PD1_int_D_predict_test,TPM_PD1_int_D_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_int_D_predict_test$xgb_TPM_PD1_int_D_predict_test,TPM_PD1_int_D_test$TPM)

# K562_PD1_int_D_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_int_D_predict_test$xgb_TPM_PD1_int_D_predict_test,TPM_PD1_int_D_test$TPM)
# 
# write.table(K562_PD1_int_D_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_int_D_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_int_D,TPM_PD1_int_D_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_neg_A ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_neg_A")
TPM_PD1_neg_A <- read.delim("~/Analysis/data/TILs/PD1_neg_A.csv",sep = ";",dec = ",")
TPM_PD1_neg_A <- merge(TPM_PD1_neg_A,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_neg_A_param <- ddply(TPM_PD1_neg_A,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_neg_A_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_neg_A$TPM <- 10^(TPM_PD1_neg_A$TPM)
TPM_PD1_neg_A <- data.frame("ID"=TPM_PD1_neg_A$ensembl_gene_id ,"TPM"=TPM_PD1_neg_A$TPM)
TPM_PD1_neg_A <- ddply(TPM_PD1_neg_A,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_neg_A <- merge(TPM_PD1_neg_A,TPM_PD1_neg_A_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_neg_A$TPM <- log10(TPM_PD1_neg_A$TPM)

TPM_PD1_neg_A[is.na(TPM_PD1_neg_A)]= 0
TPM_PD1_neg_A <- subset(TPM_PD1_neg_A, TPM_PD1_neg_A$TPM>=-1)
print(dim(TPM_PD1_neg_A)) # 13729  2682

rownames(TPM_PD1_neg_A) <- TPM_PD1_neg_A$ID
TPM_PD1_neg_A$ID <- NULL
TPM_PD1_neg_A$ensembl_gene_id <- NULL
TPM_PD1_neg_A$gene_biotype <- NULL
TPM_PD1_neg_A$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_neg_A[is.na(TPM_PD1_neg_A)]=0
TPM_PD1_neg_A <- TPM_PD1_neg_A[, -caret::nearZeroVar(TPM_PD1_neg_A, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_neg_A)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_neg_A), 0.8 * nrow(TPM_PD1_neg_A))
TPM_PD1_neg_A_train <- TPM_PD1_neg_A[train_row,]
TPM_PD1_neg_A_test <- TPM_PD1_neg_A[-train_row,]

dim(TPM_PD1_neg_A_train)
dim(TPM_PD1_neg_A_test)

write.table(TPM_PD1_neg_A_train,"~/Analysis/split_dataset/TPM_PD1_neg_A_train",sep=";")
write.table(TPM_PD1_neg_A_test,"~/Analysis/split_dataset/TPM_PD1_neg_A_test",sep=";")

# TPM_PD1_neg_A_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_A_train",sep=";")
# TPM_PD1_neg_A_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_A_test",sep=";")



## Training the XGB model ##
print("modeling PD1_neg_A")
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
xgb_TPM_PD1_neg_A <- train(TPM~.,
                      data=TPM_PD1_neg_A_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_neg_A)

saveRDS(xgb_TPM_PD1_neg_A, "~/Analysis/models/xgb_TPM_PD1_neg_A_24_06_2022.RDS")
# xgb_TPM_PD1_neg_A <- readRDS( "~/Analysis/models/xgb_TPM_PD1_neg_A_24_06_2022.RDS")

xgb_TPM_PD1_neg_A_predict_test <- predict(xgb_TPM_PD1_neg_A,TPM_PD1_neg_A_test)
xgb_TPM_PD1_neg_A_predict_test <- data.frame(xgb_TPM_PD1_neg_A_predict_test)

print(cor(xgb_TPM_PD1_neg_A_predict_test$xgb_TPM_PD1_neg_A_predict_test,TPM_PD1_neg_A_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_neg_A_predict_test$xgb_TPM_PD1_neg_A_predict_test,TPM_PD1_neg_A_test$TPM)

# K562_PD1_neg_A_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_neg_A_predict_test$xgb_TPM_PD1_neg_A_predict_test,TPM_PD1_neg_A_test$TPM)
# 
# write.table(K562_PD1_neg_A_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_neg_A_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_neg_A,TPM_PD1_neg_A_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_neg_B ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_neg_B")
TPM_PD1_neg_B <- read.delim("~/Analysis/data/TILs/PD1_neg_B.csv",sep = ";",dec = ",")
TPM_PD1_neg_B <- merge(TPM_PD1_neg_B,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_neg_B_param <- ddply(TPM_PD1_neg_B,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_neg_B_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_neg_B$TPM <- 10^(TPM_PD1_neg_B$TPM)
TPM_PD1_neg_B <- data.frame("ID"=TPM_PD1_neg_B$ensembl_gene_id ,"TPM"=TPM_PD1_neg_B$TPM)
TPM_PD1_neg_B <- ddply(TPM_PD1_neg_B,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_neg_B <- merge(TPM_PD1_neg_B,TPM_PD1_neg_B_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_neg_B$TPM <- log10(TPM_PD1_neg_B$TPM)

TPM_PD1_neg_B[is.na(TPM_PD1_neg_B)]= 0
TPM_PD1_neg_B <- subset(TPM_PD1_neg_B, TPM_PD1_neg_B$TPM>=-1)
print(dim(TPM_PD1_neg_B)) # 13729  2682

rownames(TPM_PD1_neg_B) <- TPM_PD1_neg_B$ID
TPM_PD1_neg_B$ID <- NULL
TPM_PD1_neg_B$ensembl_gene_id <- NULL
TPM_PD1_neg_B$gene_biotype <- NULL
TPM_PD1_neg_B$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_neg_B[is.na(TPM_PD1_neg_B)]=0
TPM_PD1_neg_B <- TPM_PD1_neg_B[, -caret::nearZeroVar(TPM_PD1_neg_B, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_neg_B)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_neg_B), 0.8 * nrow(TPM_PD1_neg_B))
TPM_PD1_neg_B_train <- TPM_PD1_neg_B[train_row,]
TPM_PD1_neg_B_test <- TPM_PD1_neg_B[-train_row,]

dim(TPM_PD1_neg_B_train)
dim(TPM_PD1_neg_B_test)

write.table(TPM_PD1_neg_B_train,"~/Analysis/split_dataset/TPM_PD1_neg_B_train",sep=";")
write.table(TPM_PD1_neg_B_test,"~/Analysis/split_dataset/TPM_PD1_neg_B_test",sep=";")

# TPM_PD1_neg_B_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_B_train",sep=";")
# TPM_PD1_neg_B_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_B_test",sep=";")



## Training the XGB model ##
print("modeling PD1_neg_B")
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
xgb_TPM_PD1_neg_B <- train(TPM~.,
                      data=TPM_PD1_neg_B_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_neg_B)

saveRDS(xgb_TPM_PD1_neg_B, "~/Analysis/models/xgb_TPM_PD1_neg_B_24_06_2022.RDS")
# xgb_TPM_PD1_neg_B <- readRDS( "~/Analysis/models/xgb_TPM_PD1_neg_B_24_06_2022.RDS")

xgb_TPM_PD1_neg_B_predict_test <- predict(xgb_TPM_PD1_neg_B,TPM_PD1_neg_B_test)
xgb_TPM_PD1_neg_B_predict_test <- data.frame(xgb_TPM_PD1_neg_B_predict_test)

print(cor(xgb_TPM_PD1_neg_B_predict_test$xgb_TPM_PD1_neg_B_predict_test,TPM_PD1_neg_B_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_neg_B_predict_test$xgb_TPM_PD1_neg_B_predict_test,TPM_PD1_neg_B_test$TPM)

# K562_PD1_neg_B_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_neg_B_predict_test$xgb_TPM_PD1_neg_B_predict_test,TPM_PD1_neg_B_test$TPM)
# 
# write.table(K562_PD1_neg_B_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_neg_B_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_neg_B,TPM_PD1_neg_B_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_neg_C ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_neg_C")
TPM_PD1_neg_C <- read.delim("~/Analysis/data/TILs/PD1_neg_C.csv",sep = ";",dec = ",")
TPM_PD1_neg_C <- merge(TPM_PD1_neg_C,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_neg_C_param <- ddply(TPM_PD1_neg_C,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_neg_C_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_neg_C$TPM <- 10^(TPM_PD1_neg_C$TPM)
TPM_PD1_neg_C <- data.frame("ID"=TPM_PD1_neg_C$ensembl_gene_id ,"TPM"=TPM_PD1_neg_C$TPM)
TPM_PD1_neg_C <- ddply(TPM_PD1_neg_C,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_neg_C <- merge(TPM_PD1_neg_C,TPM_PD1_neg_C_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_neg_C$TPM <- log10(TPM_PD1_neg_C$TPM)

TPM_PD1_neg_C[is.na(TPM_PD1_neg_C)]= 0
TPM_PD1_neg_C <- subset(TPM_PD1_neg_C, TPM_PD1_neg_C$TPM>=-1)
print(dim(TPM_PD1_neg_C)) # 13729  2682

rownames(TPM_PD1_neg_C) <- TPM_PD1_neg_C$ID
TPM_PD1_neg_C$ID <- NULL
TPM_PD1_neg_C$ensembl_gene_id <- NULL
TPM_PD1_neg_C$gene_biotype <- NULL
TPM_PD1_neg_C$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_neg_C[is.na(TPM_PD1_neg_C)]=0
TPM_PD1_neg_C <- TPM_PD1_neg_C[, -caret::nearZeroVar(TPM_PD1_neg_C, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_neg_C)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_neg_C), 0.8 * nrow(TPM_PD1_neg_C))
TPM_PD1_neg_C_train <- TPM_PD1_neg_C[train_row,]
TPM_PD1_neg_C_test <- TPM_PD1_neg_C[-train_row,]

dim(TPM_PD1_neg_C_train)
dim(TPM_PD1_neg_C_test)

write.table(TPM_PD1_neg_C_train,"~/Analysis/split_dataset/TPM_PD1_neg_C_train",sep=";")
write.table(TPM_PD1_neg_C_test,"~/Analysis/split_dataset/TPM_PD1_neg_C_test",sep=";")

# TPM_PD1_neg_C_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_C_train",sep=";")
# TPM_PD1_neg_C_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_C_test",sep=";")



## Training the XGB model ##
print("modeling PD1_neg_C")
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
xgb_TPM_PD1_neg_C <- train(TPM~.,
                      data=TPM_PD1_neg_C_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_neg_C)

saveRDS(xgb_TPM_PD1_neg_C, "~/Analysis/models/xgb_TPM_PD1_neg_C_24_06_2022.RDS")
# xgb_TPM_PD1_neg_C <- readRDS( "~/Analysis/models/xgb_TPM_PD1_neg_C_24_06_2022.RDS")

xgb_TPM_PD1_neg_C_predict_test <- predict(xgb_TPM_PD1_neg_C,TPM_PD1_neg_C_test)
xgb_TPM_PD1_neg_C_predict_test <- data.frame(xgb_TPM_PD1_neg_C_predict_test)

print(cor(xgb_TPM_PD1_neg_C_predict_test$xgb_TPM_PD1_neg_C_predict_test,TPM_PD1_neg_C_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_neg_C_predict_test$xgb_TPM_PD1_neg_C_predict_test,TPM_PD1_neg_C_test$TPM)

# K562_PD1_neg_C_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_neg_C_predict_test$xgb_TPM_PD1_neg_C_predict_test,TPM_PD1_neg_C_test$TPM)
# 
# write.table(K562_PD1_neg_C_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_neg_C_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_neg_C,TPM_PD1_neg_C_param)

















###__________________________________________________________________________________________###
###_________________________________________ PD1_neg_D ________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for PD1_neg_D")
TPM_PD1_neg_D <- read.delim("~/Analysis/data/TILs/PD1_neg_D.csv",sep = ";",dec = ",")
TPM_PD1_neg_D <- merge(TPM_PD1_neg_D,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_PD1_neg_D_param <- ddply(TPM_PD1_neg_D,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_PD1_neg_D_param$TPM <- NULL

print("integrating TPM per gene")
TPM_PD1_neg_D$TPM <- 10^(TPM_PD1_neg_D$TPM)
TPM_PD1_neg_D <- data.frame("ID"=TPM_PD1_neg_D$ensembl_gene_id ,"TPM"=TPM_PD1_neg_D$TPM)
TPM_PD1_neg_D <- ddply(TPM_PD1_neg_D,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_PD1_neg_D <- merge(TPM_PD1_neg_D,TPM_PD1_neg_D_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_PD1_neg_D$TPM <- log10(TPM_PD1_neg_D$TPM)

TPM_PD1_neg_D[is.na(TPM_PD1_neg_D)]= 0
TPM_PD1_neg_D <- subset(TPM_PD1_neg_D, TPM_PD1_neg_D$TPM>=-1)
print(dim(TPM_PD1_neg_D)) # 13729  2682

rownames(TPM_PD1_neg_D) <- TPM_PD1_neg_D$ID
TPM_PD1_neg_D$ID <- NULL
TPM_PD1_neg_D$ensembl_gene_id <- NULL
TPM_PD1_neg_D$gene_biotype <- NULL
TPM_PD1_neg_D$gene_name <- NULL

## Feature selection ##

print("removing useless columns")
TPM_PD1_neg_D[is.na(TPM_PD1_neg_D)]=0
TPM_PD1_neg_D <- TPM_PD1_neg_D[, -caret::nearZeroVar(TPM_PD1_neg_D, allowParallel = TRUE, uniqueCut = 0.1)]
print(dim(TPM_PD1_neg_D)) # 13729  2637

gc()


## test / train sets ##
train_row <- sample(1:nrow(TPM_PD1_neg_D), 0.8 * nrow(TPM_PD1_neg_D))
TPM_PD1_neg_D_train <- TPM_PD1_neg_D[train_row,]
TPM_PD1_neg_D_test <- TPM_PD1_neg_D[-train_row,]

dim(TPM_PD1_neg_D_train)
dim(TPM_PD1_neg_D_test)

write.table(TPM_PD1_neg_D_train,"~/Analysis/split_dataset/TPM_PD1_neg_D_train",sep=";")
write.table(TPM_PD1_neg_D_test,"~/Analysis/split_dataset/TPM_PD1_neg_D_test",sep=";")

# TPM_PD1_neg_D_train <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_D_train",sep=";")
# TPM_PD1_neg_D_test <- read.delim("~/Analysis/split_dataset/TPM_PD1_neg_D_test",sep=";")



## Training the XGB model ##
print("modeling PD1_neg_D")
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
xgb_TPM_PD1_neg_D <- train(TPM~.,
                      data=TPM_PD1_neg_D_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid,
                      na.action = na.omit,
                      nthread=80,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_PD1_neg_D)

saveRDS(xgb_TPM_PD1_neg_D, "~/Analysis/models/xgb_TPM_PD1_neg_D_24_06_2022.RDS")
# xgb_TPM_PD1_neg_D <- readRDS( "~/Analysis/models/xgb_TPM_PD1_neg_D_24_06_2022.RDS")

xgb_TPM_PD1_neg_D_predict_test <- predict(xgb_TPM_PD1_neg_D,TPM_PD1_neg_D_test)
xgb_TPM_PD1_neg_D_predict_test <- data.frame(xgb_TPM_PD1_neg_D_predict_test)

print(cor(xgb_TPM_PD1_neg_D_predict_test$xgb_TPM_PD1_neg_D_predict_test,TPM_PD1_neg_D_test$TPM, method = "pearson", use = "complete.obs")^2)

# plot(xgb_TPM_PD1_neg_D_predict_test$xgb_TPM_PD1_neg_D_predict_test,TPM_PD1_neg_D_test$TPM)

# K562_PD1_neg_D_pred_vs_test_set <- cbind.data.frame(xgb_TPM_PD1_neg_D_predict_test$xgb_TPM_PD1_neg_D_predict_test,TPM_PD1_neg_D_test$TPM)
# 
# write.table(K562_PD1_neg_D_pred_vs_test_set,"/home/nicol01b/Analysis/K562_PD1_neg_D_pred_vs_test_set.csv",sep=";",dec=",",row.names = F)

# 
# gc()
rm(TPM_PD1_neg_D,TPM_PD1_neg_D_param)

























