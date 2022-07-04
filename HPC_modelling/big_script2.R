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
#install.packages("viridis")
library(viridis)


setwd("~/Analysis/")

## Biomart ##
# ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
# tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)


## parameters ##
print("importing lib")
#Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
# dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###___________________________________________CD8_Teff_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
# print("preping data")
# TPM_CD8_Teff <- read.delim("~/Analysis/data/RNA/CD8_Teff_TPM_log10.csv",sep=";",dec=",")
# 
# TPM_CD8_Teff <- merge(TPM_CD8_Teff,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD8_Teff_param <- ddply(TPM_CD8_Teff,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD8_Teff_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD8_Teff$TPM <- 10^(TPM_CD8_Teff$TPM)
# TPM_CD8_Teff <- data.frame("ID"=TPM_CD8_Teff$ensembl_gene_id ,"TPM"=TPM_CD8_Teff$TPM)
# TPM_CD8_Teff <- ddply(TPM_CD8_Teff,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD8_Teff <- merge(TPM_CD8_Teff,TPM_CD8_Teff_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_CD8_Teff$TPM <- log10(TPM_CD8_Teff$TPM)
# 
# TPM_CD8_Teff[is.na(TPM_CD8_Teff)]= 0
# TPM_CD8_Teff <- subset(TPM_CD8_Teff, TPM_CD8_Teff$TPM>=-1)
# print(dim(TPM_CD8_Teff)) # 13729  2682
# 
# rownames(TPM_CD8_Teff) <- TPM_CD8_Teff$ID
# TPM_CD8_Teff$ID <- NULL
# TPM_CD8_Teff$ensembl_gene_id <- NULL
# TPM_CD8_Teff$gene_biotype <- NULL
# TPM_CD8_Teff$gene_name <- NULL
# 

#write.table(TPM_CD8_Teff,"~/Analysis/counts&libs/TPM_CD8_Teff_with_lib_per_gene.csv",row.names = T)
#  TPM_CD8_Teff <- read.delim("~/Analysis/counts&libs/TPM_CD8_Teff_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD8_Teff)) # 13881  7112
# 
# hist(TPM_CD8_Teff$TPM)


# test <- caret::nearZeroVar(TPM_CD8_Teff, allowParallel = TRUE, 
#                            uniqueCut = 0.1,
#                            saveMetrics = F,
#                            freqCut = dim(TPM_CD8_Teff)[1]/2)
# test$ID <- rownames(test)
# table(test$nzv)
# dim(test[test$percentUnique>0.1,])
# #test2 <- test[test$percentUnique>0.1,]
# 
# test <- data.frame("colSums"=colSums(TPM_CD8_Teff))
# test$colMeans <- colMeans(TPM_CD8_Teff)
# test$colvar <- sapply(TPM_CD8_Teff,var)
# test$sd <- sapply(TPM_CD8_Teff,sd)
# test$sd <- sapply(TPM_CD8_Teff,FUN = )
# test$ID <- rownames(test)
# 
# plot(log10(test$colMeans),log10(test$colSums))
# plot(log10(test$colvar),log10(test$colSums))
# plot(log10(test$colvar),log10(test$sd))
# 
# dim(test[test$colSums>100,])
# 
# hist((TPM_CD8_Teff$CCAATCC_UTR5))
# table(TPM_CD8_Teff$UTR5_m7G)
# 95/5
# 13881/19
# 
# 
# plot(test$nzv,log10(test$percentUnique/test$freqRatio))




## Feature selection ##
## Tn 0h ##

# print("removing useless columns")
# TPM_CD8_Teff[is.na(TPM_CD8_Teff)]=0
# #TPM_CD8_Teff_TPM <- TPM_CD8_Teff$TPM
# TPM_CD8_Teff <- TPM_CD8_Teff[, -caret::nearZeroVar(TPM_CD8_Teff, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# # col_to_keep <- data.frame("colSums"=colSums(TPM_CD8_Teff))
# # col_to_keep$ID <- rownames(col_to_keep)
# # col_to_keep <- subset(col_to_keep,col_to_keep$colSums>100)
# # 
# # dim(TPM_CD8_Teff[,col_to_keep$ID])
# # TPM_CD8_Teff <- TPM_CD8_Teff[,col_to_keep$ID]
# 
# #TPM_CD8_Teff$TPM <- TPM_CD8_Teff_TPM
# print(dim(TPM_CD8_Teff)) # 13729  2637
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(TPM_CD8_Teff), 0.8 * nrow(TPM_CD8_Teff))
# TPM_CD8_Teff_train <- TPM_CD8_Teff[train_row,]
# TPM_CD8_Teff_test <- TPM_CD8_Teff[-train_row,]
# 
# dim(TPM_CD8_Teff_train)
# dim(TPM_CD8_Teff_test)
# 
# write.table(TPM_CD8_Teff_train,"~/Analysis/TPM_CD8_Teff_train",sep=";")
# write.table(TPM_CD8_Teff_test,"TPM_CD8_Teff_test",sep=";")

TPM_CD8_Teff_train <- read.delim("~/Analysis/TPM_CD8_Teff_train",sep=";")
TPM_CD8_Teff_test <- read.delim("TPM_CD8_Teff_test",sep=";")




gc()


## Training the XGB model for Tn 0h ##
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="cv",
#                         number=10,
#                         #repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_CD8Teff <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                             max_depth = 6,
#                             colsample_bytree = 0.5,#seq(0.3, 0.5, length.out = 3),
#                             ## The values below are default values in the sklearn-api. 
#                             eta = 0.05,
#                             gamma=0,
#                             min_child_weight = 0.9,
#                             subsample = 0.5)
# 
# start_time <- Sys.time()
# xgb_TPM_CD8_Teff <- train(TPM~.,
#                       data=TPM_CD8_Teff_train,
#                       method="xgbTree",
#                       trControl=control,
#                       #metric="Rsquared",
#                       tuneGrid= xgbGrid_CD8Teff,
#                       na.action = na.omit,
#                       nthread=64,
#                       verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of 
# print(xgb_TPM_CD8_Teff)
# 
# saveRDS(xgb_TPM_CD8_Teff, "~/Analysis/models/xgb_10knround_TPM_CD8_Teff_08_10_2021.RDS")
xgb_TPM_CD8_Teff <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD8_Teff_08_10_2021.RDS")

xgb_TPM_CD8_Teff_predict_test <- predict(xgb_TPM_CD8_Teff,TPM_CD8_Teff_test)
xgb_TPM_CD8_Teff_predict_test <- data.frame(xgb_TPM_CD8_Teff_predict_test)
lm(xgb_TPM_CD8_Teff_predict_test$xgb_TPM_CD8_Teff_predict_test~TPM_CD8_Teff_test$TPM)


