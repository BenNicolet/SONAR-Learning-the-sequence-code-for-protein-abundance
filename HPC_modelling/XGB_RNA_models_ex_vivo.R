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
###___________________________________________CD8_Tn_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD8 Tn")
# TPM_CD8_Tn <- read.delim("~/Analysis/data/RNA/CD8_Tn_TPM_log10.csv",sep=";",dec=",")
# TPM_CD8_Tn <- merge(TPM_CD8_Tn,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD8_Tn_param <- ddply(TPM_CD8_Tn,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD8_Tn_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD8_Tn$TPM <- 10^(TPM_CD8_Tn$TPM)
# TPM_CD8_Tn <- data.frame("ID"=TPM_CD8_Tn$ensembl_gene_id ,"TPM"=TPM_CD8_Tn$TPM)
# TPM_CD8_Tn <- ddply(TPM_CD8_Tn,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD8_Tn <- merge(TPM_CD8_Tn,TPM_CD8_Tn_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD8_Tn$TPM <- log10(TPM_CD8_Tn$TPM)
# 
# TPM_CD8_Tn[is.na(TPM_CD8_Tn)]= 0
# TPM_CD8_Tn <- subset(TPM_CD8_Tn, TPM_CD8_Tn$TPM>=-1)
# print(dim(TPM_CD8_Tn)) # 13729  2682
# 
# rownames(TPM_CD8_Tn) <- TPM_CD8_Tn$ID
# TPM_CD8_Tn$ID <- NULL
# TPM_CD8_Tn$ensembl_gene_id <- NULL
# TPM_CD8_Tn$gene_biotype <- NULL
# TPM_CD8_Tn$gene_name <- NULL
# 
# 
# #write.table(TPM_CD8_Tn,"~/Analysis/counts&libs/TPM_CD8_Tn_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD8_Tn <- read.delim("~/Analysis/counts&libs/TPM_CD8_Tn_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD8_Tn)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD8_Tn")
# TPM_CD8_Tn[is.na(TPM_CD8_Tn)]=0
# TPM_CD8_Tn <- TPM_CD8_Tn[, -caret::nearZeroVar(TPM_CD8_Tn, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD8_Tn), 0.8 * nrow(TPM_CD8_Tn))
# TPM_CD8_Tn_train <- TPM_CD8_Tn[train_row,]
# TPM_CD8_Tn_test <- TPM_CD8_Tn[-train_row,]
# 
# dim(TPM_CD8_Tn_train)
# dim(TPM_CD8_Tn_test)
# 
# write.table(TPM_CD8_Tn_train,"~/Analysis/split_dataset/TPM_CD8_Tn_train",sep=";")
# write.table(TPM_CD8_Tn_test,"~/Analysis/split_dataset/TPM_CD8_Tn_test",sep=";")

TPM_CD8_Tn_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_train",sep=";")
TPM_CD8_Tn_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_test",sep=";")




gc()


## Training the XGB model for Tn 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8Tn <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                               max_depth = 6,
                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.05,
                               gamma=0,
                               min_child_weight = 0.9,
                               subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD8_Tn <- train(TPM~.,
                          data=TPM_CD8_Tn_train,
                          method="xgbTree",
                          trControl=control,
                          #metric="Rsquared",
                          tuneGrid= xgbGrid_CD8Tn,
                          na.action = na.omit,
                          nthread=94,
                          verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD8_Tn)

saveRDS(xgb_TPM_CD8_Tn, "~/Analysis/models/xgb_10knround_TPM_CD8_Tn_14_10_2021.RDS")
xgb_TPM_CD8_Tn <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD8_Tn_14_10_2021.RDS")

xgb_TPM_CD8_Tn_predict_test <- predict(xgb_TPM_CD8_Tn,TPM_CD8_Tn_test)
xgb_TPM_CD8_Tn_predict_test <- data.frame(xgb_TPM_CD8_Tn_predict_test)
print(lm(xgb_TPM_CD8_Tn_predict_test$xgb_TPM_CD8_Tn_predict_test~TPM_CD8_Tn_test$TPM))
plot(xgb_TPM_CD8_Tn_predict_test$xgb_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM)


print(cor(xgb_TPM_CD8_Tn_predict_test$xgb_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD8_Tn_predict_test$xgb_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM)

pred_vs_measured_CD8_Tn <- cbind.data.frame(xgb_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM)

ggplot(pred_vs_measured_CD8_Tn,aes(x=xgb_TPM_CD8_Tn_predict_test,y=TPM_CD8_Tn_test$TPM))+
  geom_pointdensity(show.legend = F, stroke=0,size=2)+
  scale_color_viridis(option = "H")+
  geom_smooth(method = "glm", color="black")+
  theme_minimal()+
  scale_x_continuous(limits = c(-1,4))+
  scale_y_continuous(limits = c(-1,4))+
  xlab("CD8 Tn predicted")+
  ylab("CD8 Tn measured")+
  theme(aspect.ratio = 1)



gc()
rm(xgb_TPM_CD8_Tn,TPM_CD8_Tn_train,TPM_CD8_Tn_test)






###__________________________________________________________________________________________###
###___________________________________________CD8_Tcm_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD8 Tcm")
# TPM_CD8_Tcm <- read.delim("~/Analysis/data/RNA/CD8_Tcm_TPM_log10.csv",sep=";",dec=",")
# TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD8_Tcm_param <- ddply(TPM_CD8_Tcm,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD8_Tcm_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD8_Tcm$TPM <- 10^(TPM_CD8_Tcm$TPM)
# TPM_CD8_Tcm <- data.frame("ID"=TPM_CD8_Tcm$ensembl_gene_id ,"TPM"=TPM_CD8_Tcm$TPM)
# TPM_CD8_Tcm <- ddply(TPM_CD8_Tcm,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,TPM_CD8_Tcm_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD8_Tcm$TPM <- log10(TPM_CD8_Tcm$TPM)
# 
# TPM_CD8_Tcm[is.na(TPM_CD8_Tcm)]= 0
# TPM_CD8_Tcm <- subset(TPM_CD8_Tcm, TPM_CD8_Tcm$TPM>=-1)
# print(dim(TPM_CD8_Tcm)) # 13729  2682
# 
# rownames(TPM_CD8_Tcm) <- TPM_CD8_Tcm$ID
# TPM_CD8_Tcm$ID <- NULL
# TPM_CD8_Tcm$ensembl_gene_id <- NULL
# TPM_CD8_Tcm$gene_biotype <- NULL
# TPM_CD8_Tcm$gene_name <- NULL
# 
# 
# #write.table(TPM_CD8_Tcm,"~/Analysis/counts&libs/TPM_CD8_Tcm_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD8_Tcm <- read.delim("~/Analysis/counts&libs/TPM_CD8_Tcm_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD8_Tcm)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD8_Tcm")
# TPM_CD8_Tcm[is.na(TPM_CD8_Tcm)]=0
# TPM_CD8_Tcm <- TPM_CD8_Tcm[, -caret::nearZeroVar(TPM_CD8_Tcm, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD8_Tcm), 0.8 * nrow(TPM_CD8_Tcm))
# TPM_CD8_Tcm_train <- TPM_CD8_Tcm[train_row,]
# TPM_CD8_Tcm_test <- TPM_CD8_Tcm[-train_row,]
# 
# dim(TPM_CD8_Tcm_train)
# dim(TPM_CD8_Tcm_test)
# 
# write.table(TPM_CD8_Tcm_train,"~/Analysis/split_dataset/TPM_CD8_Tcm_train",sep=";")
# write.table(TPM_CD8_Tcm_test,"~/Analysis/split_dataset/TPM_CD8_Tcm_test",sep=";")

TPM_CD8_Tcm_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tcm_train",sep=";")
TPM_CD8_Tcm_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tcm_test",sep=";")




gc()


## Training the XGB model for Tcm 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8Tcm <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                             max_depth = 6,
                             colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                             ## The values below are default values in the sklearn-api.
                             eta = 0.05,
                             gamma=0,
                             min_child_weight = 0.9,
                             subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD8_Tcm <- train(TPM~.,
                        data=TPM_CD8_Tcm_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD8Tcm,
                        na.action = na.omit,
                        nthread=94,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD8_Tcm)

saveRDS(xgb_TPM_CD8_Tcm, "~/Analysis/models/xgb_10knround_TPM_CD8_Tcm_14_10_2021.RDS")
xgb_TPM_CD8_Tcm <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD8_Tcm_14_10_2021.RDS")

xgb_TPM_CD8_Tcm_predict_test <- predict(xgb_TPM_CD8_Tcm,TPM_CD8_Tcm_test)
xgb_TPM_CD8_Tcm_predict_test <- data.frame(xgb_TPM_CD8_Tcm_predict_test)
print(lm(xgb_TPM_CD8_Tcm_predict_test$xgb_TPM_CD8_Tcm_predict_test~TPM_CD8_Tcm_test$TPM))


print(cor(xgb_TPM_CD8_Tcm_predict_test$xgb_TPM_CD8_Tcm_predict_test,TPM_CD8_Tcm_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD8_Tcm_predict_test$xgb_TPM_CD8_Tcm_predict_test,TPM_CD8_Tcm_test$TPM)

xgb_TPM_CD8_Tcm_imp <- data.frame(varImp(xgb_TPM_CD8_Tcm)$importance)
xgb_TPM_CD8_Tcm_imp$ID <- rownames(xgb_TPM_CD8_Tcm_imp)



gc()
rm(xgb_TPM_CD8_Tcm,TPM_CD8_Tcm_train,TPM_CD8_Tcm_test)








###__________________________________________________________________________________________###
###___________________________________________CD8_Tem_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD8 Tem")
# TPM_CD8_Tem <- read.delim("~/Analysis/data/RNA/CD8_Tem_TPM_log10.csv",sep=";",dec=",")
# TPM_CD8_Tem <- merge(TPM_CD8_Tem,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD8_Tem_param <- ddply(TPM_CD8_Tem,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD8_Tem_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD8_Tem$TPM <- 10^(TPM_CD8_Tem$TPM)
# TPM_CD8_Tem <- data.frame("ID"=TPM_CD8_Tem$ensembl_gene_id ,"TPM"=TPM_CD8_Tem$TPM)
# TPM_CD8_Tem <- ddply(TPM_CD8_Tem,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD8_Tem <- merge(TPM_CD8_Tem,TPM_CD8_Tem_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD8_Tem$TPM <- log10(TPM_CD8_Tem$TPM)
# 
# TPM_CD8_Tem[is.na(TPM_CD8_Tem)]= 0
# TPM_CD8_Tem <- subset(TPM_CD8_Tem, TPM_CD8_Tem$TPM>=-1)
# print(dim(TPM_CD8_Tem)) # 13729  2682
# 
# rownames(TPM_CD8_Tem) <- TPM_CD8_Tem$ID
# TPM_CD8_Tem$ID <- NULL
# TPM_CD8_Tem$ensembl_gene_id <- NULL
# TPM_CD8_Tem$gene_biotype <- NULL
# TPM_CD8_Tem$gene_name <- NULL
# 
# 
# #write.table(TPM_CD8_Tem,"~/Analysis/counts&libs/TPM_CD8_Tem_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD8_Tem <- read.delim("~/Analysis/counts&libs/TPM_CD8_Tem_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD8_Tem)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD8_Tem")
# TPM_CD8_Tem[is.na(TPM_CD8_Tem)]=0
# TPM_CD8_Tem <- TPM_CD8_Tem[, -caret::nearZeroVar(TPM_CD8_Tem, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD8_Tem), 0.8 * nrow(TPM_CD8_Tem))
# TPM_CD8_Tem_train <- TPM_CD8_Tem[train_row,]
# TPM_CD8_Tem_test <- TPM_CD8_Tem[-train_row,]
# 
# dim(TPM_CD8_Tem_train)
# dim(TPM_CD8_Tem_test)
# 
# write.table(TPM_CD8_Tem_train,"~/Analysis/split_dataset/TPM_CD8_Tem_train",sep=";")
# write.table(TPM_CD8_Tem_test,"~/Analysis/split_dataset/TPM_CD8_Tem_test",sep=";")

TPM_CD8_Tem_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tem_train",sep=";")
TPM_CD8_Tem_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tem_test",sep=";")




gc()


## Training the XGB model for Tem 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8Tem <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                              max_depth = 6,
                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                              ## The values below are default values in the sklearn-api.
                              eta = 0.05,
                              gamma=0,
                              min_child_weight = 0.9,
                              subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD8_Tem <- train(TPM~.,
                         data=TPM_CD8_Tem_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_CD8Tem,
                         na.action = na.omit,
                         nthread=94,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD8_Tem)

saveRDS(xgb_TPM_CD8_Tem, "~/Analysis/models/xgb_10knround_TPM_CD8_Tem_14_10_2021.RDS")
xgb_TPM_CD8_Tem <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD8_Tem_14_10_2021.RDS")

xgb_TPM_CD8_Tem_predict_test <- predict(xgb_TPM_CD8_Tem,TPM_CD8_Tem_test)
xgb_TPM_CD8_Tem_predict_test <- data.frame(xgb_TPM_CD8_Tem_predict_test)
print(lm(xgb_TPM_CD8_Tem_predict_test$xgb_TPM_CD8_Tem_predict_test~TPM_CD8_Tem_test$TPM))

print(cor(xgb_TPM_CD8_Tem_predict_test$xgb_TPM_CD8_Tem_predict_test,TPM_CD8_Tem_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD8_Tem_predict_test$xgb_TPM_CD8_Tem_predict_test,TPM_CD8_Tem_test$TPM)


gc()
rm(xgb_TPM_CD8_Tem,TPM_CD8_Tem_train,TPM_CD8_Tem_test)












###__________________________________________________________________________________________###
###___________________________________________CD4_Tn_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD4 Tn")
# TPM_CD4_Tn <- read.delim("~/Analysis/data/RNA/CD4_Tn_TPM_log10.csv",sep=";",dec=",")
# TPM_CD4_Tn <- merge(TPM_CD4_Tn,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD4_Tn_param <- ddply(TPM_CD4_Tn,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD4_Tn_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD4_Tn$TPM <- 10^(TPM_CD4_Tn$TPM)
# TPM_CD4_Tn <- data.frame("ID"=TPM_CD4_Tn$ensembl_gene_id ,"TPM"=TPM_CD4_Tn$TPM)
# TPM_CD4_Tn <- ddply(TPM_CD4_Tn,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD4_Tn <- merge(TPM_CD4_Tn,TPM_CD4_Tn_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD4_Tn$TPM <- log10(TPM_CD4_Tn$TPM)
# 
# TPM_CD4_Tn[is.na(TPM_CD4_Tn)]= 0
# TPM_CD4_Tn <- subset(TPM_CD4_Tn, TPM_CD4_Tn$TPM>=-1)
# print(dim(TPM_CD4_Tn)) # 13729  2682
# 
# rownames(TPM_CD4_Tn) <- TPM_CD4_Tn$ID
# TPM_CD4_Tn$ID <- NULL
# TPM_CD4_Tn$ensembl_gene_id <- NULL
# TPM_CD4_Tn$gene_biotype <- NULL
# TPM_CD4_Tn$gene_name <- NULL
# 
# 
# #write.table(TPM_CD4_Tn,"~/Analysis/counts&libs/TPM_CD4_Tn_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD4_Tn <- read.delim("~/Analysis/counts&libs/TPM_CD4_Tn_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD4_Tn)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD4_Tn")
# TPM_CD4_Tn[is.na(TPM_CD4_Tn)]=0
# TPM_CD4_Tn <- TPM_CD4_Tn[, -caret::nearZeroVar(TPM_CD4_Tn, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD4_Tn), 0.8 * nrow(TPM_CD4_Tn))
# TPM_CD4_Tn_train <- TPM_CD4_Tn[train_row,]
# TPM_CD4_Tn_test <- TPM_CD4_Tn[-train_row,]
# 
# dim(TPM_CD4_Tn_train)
# dim(TPM_CD4_Tn_test)
# 
# write.table(TPM_CD4_Tn_train,"~/Analysis/split_dataset/TPM_CD4_Tn_train",sep=";")
# write.table(TPM_CD4_Tn_test,"~/Analysis/split_dataset/TPM_CD4_Tn_test",sep=";")

TPM_CD4_Tn_train <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tn_train",sep=";")
TPM_CD4_Tn_test <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tn_test",sep=";")




gc()


## Training the XGB model for Tn 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4Tn <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                              max_depth = 6,
                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                              ## The values below are default values in the sklearn-api.
                              eta = 0.05,
                              gamma=0,
                              min_child_weight = 0.9,
                              subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD4_Tn <- train(TPM~.,
                         data=TPM_CD4_Tn_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_CD4Tn,
                         na.action = na.omit,
                         nthread=94,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD4_Tn)

saveRDS(xgb_TPM_CD4_Tn, "~/Analysis/models/xgb_10knround_TPM_CD4_Tn_14_10_2021.RDS")
xgb_TPM_CD4_Tn <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD4_Tn_14_10_2021.RDS")

xgb_TPM_CD4_Tn_predict_test <- predict(xgb_TPM_CD4_Tn,TPM_CD4_Tn_test)
xgb_TPM_CD4_Tn_predict_test <- data.frame(xgb_TPM_CD4_Tn_predict_test)
print(lm(xgb_TPM_CD4_Tn_predict_test$xgb_TPM_CD4_Tn_predict_test~TPM_CD4_Tn_test$TPM))

print(cor(xgb_TPM_CD4_Tn_predict_test$xgb_TPM_CD4_Tn_predict_test,TPM_CD4_Tn_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD4_Tn_predict_test$xgb_TPM_CD4_Tn_predict_test,TPM_CD4_Tn_test$TPM)


gc()
rm(xgb_TPM_CD4_Tn,TPM_CD4_Tn_train,TPM_CD4_Tn_test)









###__________________________________________________________________________________________###
###___________________________________________CD4_Tcm_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD4 Tcm")
# TPM_CD4_Tcm <- read.delim("~/Analysis/data/RNA/CD4_Tcm_TPM_log10.csv",sep=";",dec=",")
# TPM_CD4_Tcm <- merge(TPM_CD4_Tcm,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD4_Tcm_param <- ddply(TPM_CD4_Tcm,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD4_Tcm_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD4_Tcm$TPM <- 10^(TPM_CD4_Tcm$TPM)
# TPM_CD4_Tcm <- data.frame("ID"=TPM_CD4_Tcm$ensembl_gene_id ,"TPM"=TPM_CD4_Tcm$TPM)
# TPM_CD4_Tcm <- ddply(TPM_CD4_Tcm,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD4_Tcm <- merge(TPM_CD4_Tcm,TPM_CD4_Tcm_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD4_Tcm$TPM <- log10(TPM_CD4_Tcm$TPM)
# 
# TPM_CD4_Tcm[is.na(TPM_CD4_Tcm)]= 0
# TPM_CD4_Tcm <- subset(TPM_CD4_Tcm, TPM_CD4_Tcm$TPM>=-1)
# print(dim(TPM_CD4_Tcm)) # 13729  2682
# 
# rownames(TPM_CD4_Tcm) <- TPM_CD4_Tcm$ID
# TPM_CD4_Tcm$ID <- NULL
# TPM_CD4_Tcm$ensembl_gene_id <- NULL
# TPM_CD4_Tcm$gene_biotype <- NULL
# TPM_CD4_Tcm$gene_name <- NULL
# 
# 
# #write.table(TPM_CD4_Tcm,"~/Analysis/counts&libs/TPM_CD4_Tcm_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD4_Tcm <- read.delim("~/Analysis/counts&libs/TPM_CD4_Tcm_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD4_Tcm)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD4_Tcm")
# TPM_CD4_Tcm[is.na(TPM_CD4_Tcm)]=0
# TPM_CD4_Tcm <- TPM_CD4_Tcm[, -caret::nearZeroVar(TPM_CD4_Tcm, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD4_Tcm), 0.8 * nrow(TPM_CD4_Tcm))
# TPM_CD4_Tcm_train <- TPM_CD4_Tcm[train_row,]
# TPM_CD4_Tcm_test <- TPM_CD4_Tcm[-train_row,]
# 
# dim(TPM_CD4_Tcm_train)
# dim(TPM_CD4_Tcm_test)
# 
# write.table(TPM_CD4_Tcm_train,"~/Analysis/split_dataset/TPM_CD4_Tcm_train",sep=";")
# write.table(TPM_CD4_Tcm_test,"~/Analysis/split_dataset/TPM_CD4_Tcm_test",sep=";")

TPM_CD4_Tcm_train <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tcm_train",sep=";")
TPM_CD4_Tcm_test <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tcm_test",sep=";")




gc()


## Training the XGB model for Tcm 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4Tcm <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                             max_depth = 6,
                             colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                             ## The values below are default values in the sklearn-api.
                             eta = 0.05,
                             gamma=0,
                             min_child_weight = 0.9,
                             subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD4_Tcm <- train(TPM~.,
                        data=TPM_CD4_Tcm_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD4Tcm,
                        na.action = na.omit,
                        nthread=94,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD4_Tcm)

saveRDS(xgb_TPM_CD4_Tcm, "~/Analysis/models/xgb_10knround_TPM_CD4_Tcm_14_10_2021.RDS")
xgb_TPM_CD4_Tcm <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD4_Tcm_14_10_2021.RDS")

xgb_TPM_CD4_Tcm_predict_test <- predict(xgb_TPM_CD4_Tcm,TPM_CD4_Tcm_test)
xgb_TPM_CD4_Tcm_predict_test <- data.frame(xgb_TPM_CD4_Tcm_predict_test)
print(lm(xgb_TPM_CD4_Tcm_predict_test$xgb_TPM_CD4_Tcm_predict_test~TPM_CD4_Tcm_test$TPM))

print(cor(xgb_TPM_CD4_Tcm_predict_test$xgb_TPM_CD4_Tcm_predict_test,TPM_CD4_Tcm_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD4_Tcm_predict_test$xgb_TPM_CD4_Tcm_predict_test,TPM_CD4_Tcm_test$TPM)


gc()
rm(xgb_TPM_CD4_Tcm,TPM_CD4_Tcm_train,TPM_CD4_Tcm_test)










###__________________________________________________________________________________________###
###___________________________________________CD4_Tem_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD4 Tem")
# TPM_CD4_Tem <- read.delim("~/Analysis/data/RNA/CD4_Tem_TPM_log10.csv",sep=";",dec=",")
# TPM_CD4_Tem <- merge(TPM_CD4_Tem,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_CD4_Tem_param <- ddply(TPM_CD4_Tem,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_CD4_Tem_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_CD4_Tem$TPM <- 10^(TPM_CD4_Tem$TPM)
# TPM_CD4_Tem <- data.frame("ID"=TPM_CD4_Tem$ensembl_gene_id ,"TPM"=TPM_CD4_Tem$TPM)
# TPM_CD4_Tem <- ddply(TPM_CD4_Tem,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD4_Tem <- merge(TPM_CD4_Tem,TPM_CD4_Tem_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# TPM_CD4_Tem$TPM <- log10(TPM_CD4_Tem$TPM)
# 
# TPM_CD4_Tem[is.na(TPM_CD4_Tem)]= 0
# TPM_CD4_Tem <- subset(TPM_CD4_Tem, TPM_CD4_Tem$TPM>=-1)
# print(dim(TPM_CD4_Tem)) # 13729  2682
# 
# rownames(TPM_CD4_Tem) <- TPM_CD4_Tem$ID
# TPM_CD4_Tem$ID <- NULL
# TPM_CD4_Tem$ensembl_gene_id <- NULL
# TPM_CD4_Tem$gene_biotype <- NULL
# TPM_CD4_Tem$gene_name <- NULL
# 
# 
# #write.table(TPM_CD4_Tem,"~/Analysis/counts&libs/TPM_CD4_Tem_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD4_Tem <- read.delim("~/Analysis/counts&libs/TPM_CD4_Tem_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD4_Tem)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD4_Tem")
# TPM_CD4_Tem[is.na(TPM_CD4_Tem)]=0
# TPM_CD4_Tem <- TPM_CD4_Tem[, -caret::nearZeroVar(TPM_CD4_Tem, allowParallel = TRUE, uniqueCut = 0.1)]
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD4_Tem), 0.8 * nrow(TPM_CD4_Tem))
# TPM_CD4_Tem_train <- TPM_CD4_Tem[train_row,]
# TPM_CD4_Tem_test <- TPM_CD4_Tem[-train_row,]
# 
# dim(TPM_CD4_Tem_train)
# dim(TPM_CD4_Tem_test)
# 
# write.table(TPM_CD4_Tem_train,"~/Analysis/split_dataset/TPM_CD4_Tem_train",sep=";")
# write.table(TPM_CD4_Tem_test,"~/Analysis/split_dataset/TPM_CD4_Tem_test",sep=";")

TPM_CD4_Tem_train <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tem_train",sep=";")
TPM_CD4_Tem_test <- read.delim("~/Analysis/split_dataset/TPM_CD4_Tem_test",sep=";")




gc()


## Training the XGB model for Tem 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD4Tem <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                             max_depth = 6,
                             colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                             ## The values below are default values in the sklearn-api.
                             eta = 0.05,
                             gamma=0,
                             min_child_weight = 0.9,
                             subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD4_Tem <- train(TPM~.,
                        data=TPM_CD4_Tem_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD4Tem,
                        na.action = na.omit,
                        nthread=94,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD4_Tem)

saveRDS(xgb_TPM_CD4_Tem, "~/Analysis/models/xgb_10knround_TPM_CD4_Tem_14_10_2021.RDS")
xgb_TPM_CD4_Tem <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD4_Tem_14_10_2021.RDS")

xgb_TPM_CD4_Tem_predict_test <- predict(xgb_TPM_CD4_Tem,TPM_CD4_Tem_test)
xgb_TPM_CD4_Tem_predict_test <- data.frame(xgb_TPM_CD4_Tem_predict_test)
print(lm(xgb_TPM_CD4_Tem_predict_test$xgb_TPM_CD4_Tem_predict_test~TPM_CD4_Tem_test$TPM))

print(cor(xgb_TPM_CD4_Tem_predict_test$xgb_TPM_CD4_Tem_predict_test,TPM_CD4_Tem_test$TPM, method = "pearson", use = "complete.obs")^2)
plot(xgb_TPM_CD4_Tem_predict_test$xgb_TPM_CD4_Tem_predict_test,TPM_CD4_Tem_test$TPM)


gc()
rm(xgb_TPM_CD4_Tem,TPM_CD4_Tem_train,TPM_CD4_Tem_test)


















###__________________________________________________________________________________________###
###____________________________________T cell activation_____________________________________###
###__________________________________________________________________________________________###

#ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# attributes_ens <- listAttributes(ensembl)
# View(attributes_ens)
#IDs_genenames_coding <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","uniprotswissprot","gene_biotype"), mart = ensembl, filters = "with_ccds", values = TRUE)
#IDs_genenames_coding <- IDs_genenames_coding[!duplicated(IDs_genenames_coding$uniprotswissprot),]


## parameters ##
print("importing library")
# Sequence_parameters_Geiger <- read.delim("libraries/Sequence_parameters_for_Geiger_data_prot_libv2_09092021.csv",sep = ";",dec=",")




###__________________________________________________________________________________________###
###___________________________________________CD4_T0h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
# print("importing and preparing data")
# CN_0h <- read.delim("data/protein/CD4_Tn_0h_CN_log10.csv",sep = ";",dec = ",")
# CN_0h <- merge(CN_0h,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(CN_0h)) # 6993 6803
# rownames(CN_0h) <- CN_0h$ID
# CN_0h$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# CN_0h[is.na(CN_0h)]=0
# CN_0h <- CN_0h[, -nearZeroVar(CN_0h, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_0h)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(CN_0h), 0.8 * nrow(CN_0h))
# CN_0h_train <- CN_0h[train_row,]
# CN_0h_test <- CN_0h[-train_row,]
# 
# print(dim(CN_0h_train))
# print(dim(CN_0h_test))
# 
# write.table(CN_0h_train,"~/Analysis/split_dataset/CN_0h_train.csv",sep=";")
# write.table(CN_0h_test,"~/Analysis/split_dataset/CN_0h_test.csv",sep=";")

CN_0h_train <- read.delim("~/Analysis/split_dataset/CN_0h_train.csv",sep=";")
CN_0h_test <- read.delim("~/Analysis/split_dataset/CN_0h_test.csv",sep=";")




gc()


## Training the RF model for Tn 0h ##
print("Modeling Tn 0h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid_CD4_0h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                              max_depth = 6,
                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                              ## The values below are default values in the sklearn-api.
                              eta = 0.05,
                              gamma=0,
                              min_child_weight = 0.9,
                              subsample = 1)


start_time <- Sys.time()
xgb_CN_Tn_0h <- train(CN~.,
                      data=CN_0h_train,
                      method="xgbTree",
                      trControl=control,
                      tuneGrid= xgbGrid_CD4_0h,
                      na.action = na.omit,
                      nthread=94,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) 
print(xgb_CN_Tn_0h)


saveRDS(xgb_CN_Tn_0h, "models/protein/xgb_CN_Tn_0h_10k_15_10_2021.RDS")
xgb_CN_Tn_0h <- readRDS("~/Analysis/models/protein/xgb_CN_Tn_0h_10k_15_10_2021.RDS")

xgb_CN_Tn_0h_predict_test <- predict(xgb_CN_Tn_0h,CN_0h_test)
xgb_CN_Tn_0h_predict_test <- data.frame(xgb_CN_Tn_0h_predict_test)
print(lm(xgb_CN_Tn_0h_predict_test$xgb_CN_Tn_0h_predict_test~CN_0h_test$CN))

print(cor(xgb_CN_Tn_0h_predict_test$xgb_CN_Tn_0h_predict_test,CN_0h_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CN_Tn_0h_predict_test$xgb_CN_Tn_0h_predict_test,CN_0h_test$CN)



Var_imp_xgb_CN_Tn_0h <- data.frame(varImp(xgb_CN_Tn_0h)$importance)
Var_imp_xgb_CN_Tn_0h$ID <- rownames(Var_imp_xgb_CN_Tn_0h)

gc()
rm(xgb_CN_Tn_0h,CN_0h_train,CN_0h_test)



###__________________________________________________________________________________________###
###___________________________________________CD4_T6h________________________________________###
###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# print("importing and preparing data CN Tn 6h")
# CN_6h <- read.delim("data/protein/CD4_Tn_6h_CN_log10.csv",sep = ";",dec = ",")
# CN_6h <- merge(CN_6h,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(CN_6h)) # 
# rownames(CN_6h) <- CN_6h$ID
# CN_6h$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# CN_6h[is.na(CN_6h)]=0
# CN_6h <- CN_6h[, -nearZeroVar(CN_6h, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_6h)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(CN_6h), 0.8 * nrow(CN_6h))
# CN_6h_train <- CN_6h[train_row,]
# CN_6h_test <- CN_6h[-train_row,]
# 
# print(dim(CN_6h_train))
# print(dim(CN_6h_test))
# 
# write.table(CN_6h_train,"~/Analysis/split_dataset/CN_6h_train.csv",sep=";")
# write.table(CN_6h_test,"~/Analysis/split_dataset/CN_6h_test.csv",sep=";")

CN_6h_train <- read.delim("~/Analysis/split_dataset/CN_6h_train.csv",sep=";")
CN_6h_test <- read.delim("~/Analysis/split_dataset/CN_6h_test.csv",sep=";")




gc()


## Training the RF model for Tn 6h ##
print("Modeling Tn 6h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid_CD4_6h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                              max_depth = 6,
                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                              ## The values below are default values in the sklearn-api.
                              eta = 0.05,
                              gamma=0,
                              min_child_weight = 0.9,
                              subsample = 1)


start_time <- Sys.time()
xgb_CN_Tn_6h <- train(CN~.,
                      data=CN_6h_train,
                      method="xgbTree",
                      trControl=control,
                      tuneGrid= xgbGrid_CD4_6h,
                      na.action = na.omit,
                      nthread=94,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(xgb_CN_Tn_6h)


saveRDS(xgb_CN_Tn_6h, "models/protein/xgb_CN_Tn_6h_10k_15_10_2021.RDS")
xgb_CN_Tn_6h <- readRDS("~/Analysis/models/protein/xgb_CN_Tn_6h_10k_15_10_2021.RDS")


xgb_CN_Tn_6h_predict_test <- predict(xgb_CN_Tn_6h,CN_6h_test)
xgb_CN_Tn_6h_predict_test <- data.frame(xgb_CN_Tn_6h_predict_test)
print(lm(xgb_CN_Tn_6h_predict_test$xgb_CN_Tn_6h_predict_test~CN_6h_test$CN))

print(cor(xgb_CN_Tn_6h_predict_test$xgb_CN_Tn_6h_predict_test,CN_6h_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CN_Tn_6h_predict_test$xgb_CN_Tn_6h_predict_test,CN_6h_test$CN)


Var_imp_xgb_CN_Tn_6h <- data.frame(varImp(xgb_CN_Tn_6h)$importance)
Var_imp_xgb_CN_Tn_6h$ID <- rownames(Var_imp_xgb_CN_Tn_6h)

gc()
rm(xgb_CN_Tn_6h,CN_6h_train,CN_6h_test)





###__________________________________________________________________________________________###
###___________________________________________CD4_T12h________________________________________###
###__________________________________________________________________________________________###

# 
# ## prep data ##
# print("importing and preparing data CN Tn 12h")
# CN_12h <- read.delim("data/protein/CD4_Tn_12h_CN_log10.csv",sep = ";",dec = ",")
# CN_12h <- merge(CN_12h,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(CN_12h)) # 
# rownames(CN_12h) <- CN_12h$ID
# CN_12h$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# CN_12h[is.na(CN_12h)]=0
# CN_12h <- CN_12h[, -nearZeroVar(CN_12h, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_12h)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(CN_12h), 0.8 * nrow(CN_12h))
# CN_12h_train <- CN_12h[train_row,]
# CN_12h_test <- CN_12h[-train_row,]
# 
# print(dim(CN_12h_train))
# print(dim(CN_12h_test))
# 
# write.table(CN_12h_train,"~/Analysis/split_dataset/CN_12h_train.csv",sep=";")
# write.table(CN_12h_test,"~/Analysis/split_dataset/CN_12h_test.csv",sep=";")

CN_12h_train <- read.delim("~/Analysis/split_dataset/CN_12h_train.csv",sep=";")
CN_12h_test <- read.delim("~/Analysis/split_dataset/CN_12h_test.csv",sep=";")




gc()


## Training the RF model for Tn 12h ##
print("Modeling Tn 12h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid_CD4_12h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                               max_depth = 6,
                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.05,
                               gamma=0,
                               min_child_weight = 0.9,
                               subsample = 1)


start_time <- Sys.time()
xgb_CN_Tn_12h <- train(CN~.,
                       data=CN_12h_train,
                       method="xgbTree",
                       trControl=control,
                       tuneGrid= xgbGrid_CD4_12h,
                       na.action = na.omit,
                       nthread=94,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(xgb_CN_Tn_12h)


saveRDS(xgb_CN_Tn_12h, "models/protein/xgb_CN_Tn_12h_10k_15_10_2021.RDS")
xgb_CN_Tn_12h <- readRDS("~/Analysis/models/protein/xgb_CN_Tn_12h_10k_15_10_2021.RDS")

xgb_CN_Tn_12h_predict_test <- predict(xgb_CN_Tn_12h,CN_12h_test)
xgb_CN_Tn_12h_predict_test <- data.frame(xgb_CN_Tn_12h_predict_test)
print(lm(xgb_CN_Tn_12h_predict_test$xgb_CN_Tn_12h_predict_test~CN_12h_test$CN))

print(cor(xgb_CN_Tn_12h_predict_test$xgb_CN_Tn_12h_predict_test,CN_12h_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CN_Tn_12h_predict_test$xgb_CN_Tn_12h_predict_test,CN_12h_test$CN)

gc()
rm(xgb_CN_Tn_12h,CN_12h_train,CN_12h_test)




###__________________________________________________________________________________________###
###___________________________________________CD4_T24h________________________________________###
###__________________________________________________________________________________________###

# 
# ## prep data ##
# print("importing and preparing data CN Tn 24h")
# CN_24h <- read.delim("data/protein/CD4_Tn_24h_CN_log10.csv",sep = ";",dec = ",")
# CN_24h <- merge(CN_24h,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(CN_24h)) # 
# rownames(CN_24h) <- CN_24h$ID
# CN_24h$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# CN_24h[is.na(CN_24h)]=0
# CN_24h <- CN_24h[, -nearZeroVar(CN_24h, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_24h)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(CN_24h), 0.8 * nrow(CN_24h))
# CN_24h_train <- CN_24h[train_row,]
# CN_24h_test <- CN_24h[-train_row,]
# 
# print(dim(CN_24h_train))
# print(dim(CN_24h_test))
# 
# write.table(CN_24h_train,"~/Analysis/split_dataset/CN_24h_train.csv",sep=";")
# write.table(CN_24h_test,"~/Analysis/split_dataset/CN_24h_test.csv",sep=";")

CN_24h_train <- read.delim("~/Analysis/split_dataset/CN_24h_train.csv",sep=";")
CN_24h_test <- read.delim("~/Analysis/split_dataset/CN_24h_test.csv",sep=";")




gc()


## Training the RF model for Tn 24h ##
print("Modeling Tn 24h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid_CD4_24h <- expand.grid(nrounds = 1000,  # this is n_estimators in the python code above
                               max_depth = 6,
                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.05,
                               gamma=0,
                               min_child_weight = 0.9,
                               subsample = 1)


start_time <- Sys.time()
xgb_CN_Tn_24h <- train(CN~.,
                       data=CN_24h_train,
                       method="xgbTree",
                       trControl=control,
                       tuneGrid= xgbGrid_CD4_24h,
                       na.action = na.omit,
                       nthread=94,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(xgb_CN_Tn_24h)


saveRDS(xgb_CN_Tn_24h, "models/protein/xgb_CN_Tn_24h_10k_15_10_2021.RDS")
xgb_CN_Tn_24h <- readRDS("~/Analysis/models/protein/xgb_CN_Tn_24h_10k_15_10_2021.RDS")

xgb_CN_Tn_24h_predict_test <- predict(xgb_CN_Tn_24h,CN_24h_test)
xgb_CN_Tn_24h_predict_test <- data.frame(xgb_CN_Tn_24h_predict_test)
print(lm(xgb_CN_Tn_24h_predict_test$xgb_CN_Tn_24h_predict_test~CN_24h_test$CN))

print(cor(xgb_CN_Tn_24h_predict_test$xgb_CN_Tn_24h_predict_test,CN_24h_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CN_Tn_24h_predict_test$xgb_CN_Tn_24h_predict_test,CN_24h_test$CN)

gc()
rm(xgb_CN_Tn_24h,CN_24h_train,CN_24h_test)





###__________________________________________________________________________________________###
###___________________________________________CD4_T48h________________________________________###
###__________________________________________________________________________________________###

# 
# ## prep data ##
# print("importing and preparing data CN Tn 48h")
# CN_48h <- read.delim("data/protein/CD4_Tn_48h_CN_log10.csv",sep = ";",dec = ",")
# CN_48h <- merge(CN_48h,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(CN_48h)) # 
# rownames(CN_48h) <- CN_48h$ID
# CN_48h$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# CN_48h[is.na(CN_48h)]=0
# CN_48h <- CN_48h[, -nearZeroVar(CN_48h, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_48h)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(CN_48h), 0.8 * nrow(CN_48h))
# CN_48h_train <- CN_48h[train_row,]
# CN_48h_test <- CN_48h[-train_row,]
# 
# print(dim(CN_48h_train))
# print(dim(CN_48h_test))
# 
# write.table(CN_48h_train,"~/Analysis/split_dataset/CN_48h_train.csv",sep=";")
# write.table(CN_48h_test,"~/Analysis/split_dataset/CN_48h_test.csv",sep=";")

CN_48h_train <- read.delim("~/Analysis/split_dataset/CN_48h_train.csv",sep=";")
CN_48h_test <- read.delim("~/Analysis/split_dataset/CN_48h_test.csv",sep=";")




gc()


## Training the RF model for Tn 48h ##
print("Modeling Tn 48h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid_CD4_48h <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                               max_depth = 6,
                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.05,
                               gamma=0,
                               min_child_weight = 0.9,
                               subsample = 1)


start_time <- Sys.time()
xgb_CN_Tn_48h <- train(CN~.,
                       data=CN_48h_train,
                       method="xgbTree",
                       trControl=control,
                       tuneGrid= xgbGrid_CD4_48h,
                       na.action = na.omit,
                       nthread=94,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(xgb_CN_Tn_48h)


saveRDS(xgb_CN_Tn_48h, "models/protein/xgb_CN_Tn_48h_10k_15_10_2021.RDS")
xgb_CN_Tn_48h <- readRDS("~/Analysis/models/protein/xgb_CN_Tn_48h_10k_15_10_2021.RDS")


xgb_CN_Tn_48h_predict_test <- predict(xgb_CN_Tn_48h,CN_48h_test)
xgb_CN_Tn_48h_predict_test <- data.frame(xgb_CN_Tn_48h_predict_test)
print(lm(xgb_CN_Tn_48h_predict_test$xgb_CN_Tn_48h_predict_test~CN_48h_test$CN))

print(cor(xgb_CN_Tn_48h_predict_test$xgb_CN_Tn_48h_predict_test,CN_48h_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CN_Tn_48h_predict_test$xgb_CN_Tn_48h_predict_test,CN_48h_test$CN)



gc()
rm(xgb_CN_Tn_48h,CN_48h_train,CN_48h_test)











###__________________________________________________________________________________________###
###____________________________________________Ex vivo_______________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")

## Biomart ##
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)


## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###___________________________________________CD8_Teff_______________________________________###
###__________________________________________________________________________________________###

# ## prep data ##
# print("preping data for CD8 Teff")
# TPM_CD8_Teff <- read.delim("~/Analysis/data/RNA/CD8_Teff_TPM_log10.csv",sep=";",dec=",")
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
# 
# #write.table(TPM_CD8_Teff,"~/Analysis/counts&libs/TPM_CD8_Teff_with_lib_per_gene.csv",row.names = T)
# #  TPM_CD8_Teff <- read.delim("~/Analysis/counts&libs/TPM_CD8_Teff_with_lib_per_gene.csv",sep = " ")
# print(dim(TPM_CD8_Teff)) # 13881  7112
# 
# 
# ## Feature selection ##
# print("removing useless columns for CD8_Teff")
# TPM_CD8_Teff[is.na(TPM_CD8_Teff)]=0
# TPM_CD8_Teff <- TPM_CD8_Teff[, -caret::nearZeroVar(TPM_CD8_Teff, allowParallel = TRUE, uniqueCut = 0.1)]
#  
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_CD8_Teff), 0.8 * nrow(TPM_CD8_Teff))
# TPM_CD8_Teff_train <- TPM_CD8_Teff[train_row,]
# TPM_CD8_Teff_test <- TPM_CD8_Teff[-train_row,]
# 
# dim(TPM_CD8_Teff_train)
# dim(TPM_CD8_Teff_test)
# 
# write.table(TPM_CD8_Teff_train,"~/Analysis/split_dataset/TPM_CD8_Teff_train",sep=";")
# write.table(TPM_CD8_Teff_test,"~/Analysis/split_dataset/TPM_CD8_Teff_test",sep=";")

TPM_CD8_Teff_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Teff_train",sep=";")
TPM_CD8_Teff_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Teff_test",sep=";")




gc()


## Training the XGB model for Tn 0h ##
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8Teff <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                               max_depth = 6,
                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.05,
                               gamma=0,
                               min_child_weight = 0.9,
                               subsample = 1)

start_time <- Sys.time()
xgb_TPM_CD8_Teff <- train(TPM~.,
                          data=TPM_CD8_Teff_train,
                          method="xgbTree",
                          trControl=control,
                          #metric="Rsquared",
                          tuneGrid= xgbGrid_CD8Teff,
                          na.action = na.omit,
                          nthread=94,
                          verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(xgb_TPM_CD8_Teff)

saveRDS(xgb_TPM_CD8_Teff, "~/Analysis/models/xgb_10knround_TPM_CD8_Teff_10k_15_10_2021.RDS")
# xgb_TPM_CD8_Teff <- readRDS("~/Analysis/models/xgb_10knround_TPM_CD8_Teff_10k_15_10_2021.RDS")

xgb_TPM_CD8_Teff_predict_test <- predict(xgb_TPM_CD8_Teff,TPM_CD8_Teff_test)
xgb_TPM_CD8_Teff_predict_test <- data.frame(xgb_TPM_CD8_Teff_predict_test)
print(lm(xgb_TPM_CD8_Teff_predict_test$xgb_TPM_CD8_Teff_predict_test~TPM_CD8_Teff_test$TPM))

gc()
rm(xgb_TPM_CD8_Teff,TPM_CD8_Teff_train,TPM_CD8_Teff_test)




