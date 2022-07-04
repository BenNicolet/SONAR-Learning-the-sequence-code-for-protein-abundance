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
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
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
###____________________________________________CD8_Tn________________________________________###
###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
# CN_CD8_Tn <- read.delim("~/Analysis/data/protein/CD8_Tn_CN_log10.csv",sep=";",dec = ",")
# CN_CD8_Tn <- subset(CN_CD8_Tn,CN_CD8_Tn$CN>0)
# dim(CN_CD8_Tn)
# 
# 
# TPM_CD8_Tn <- read.delim("~/Analysis/data/RNA/CD8_Tn_TPM_log10.csv",sep=";",dec=",")
# TPM_CD8_Tn$TPM <- 10^TPM_CD8_Tn$TPM
# TPM_CD8_Tn <- merge(TPM_CD8_Tn,tx2gene, by.x="ID", by.y="ensembl_transcript_id")
# 
# 
# registerDoMC(10)
# 
# TPM_CD8_Tn <- merge(TPM_CD8_Tn,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
# TPM_CD8_Tn <- ddply(TPM_CD8_Tn,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
# TPM_CD8_Tn$TPM <- log10(TPM_CD8_Tn$TPM)
# #dim(TPM_CD8_Tn[duplicated(TPM_CD8_Tn$uniprotswissprot),])
# 
# 
# 
# CN_CD8_Tn <- merge(CN_CD8_Tn,TPM_CD8_Tn,by.x="ID",by.y="uniprotswissprot")
# CN_CD8_Tn$ensembl_gene_id <- NULL
# CN_CD8_Tn <- merge(CN_CD8_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_CD8_Tn[is.na(CN_CD8_Tn)]= 0
# print(dim(CN_CD8_Tn)) # 7984 7126
# 
# CN_CD8_Tn_dups <- CN_CD8_Tn[duplicated(CN_CD8_Tn$ID),]
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
# # write.table(CN_CD8_Tn_train,"~/Analysis/split_dataset/CN_CD8_Tn_with_RNA_train",sep=";")
# # write.table(CN_CD8_Tn_test,"~/Analysis/split_dataset/CN_CD8_Tn_with_RNA_test",sep=";")
# 
# CN_CD8_Tn_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tn_with_RNA_train",sep=";")
# CN_CD8_Tn_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tn_with_RNA_test",sep=";")
# # 
# # 
# # 
# # ## Training the xgb model ##
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
#                         nthread=80,
#                         verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_CD8_TN_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tn_with_RNA_18_10_2021.RDS")
# # xgb_CD8_TN_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tn_with_RNA_18_10_2021.RDS")
# # 
# # print(xgb_CD8_TN_CN)
# 
# xgb_CD8_TN_CN_predict_test <- predict(xgb_CD8_TN_CN,CN_CD8_Tn_test)
# xgb_CD8_TN_CN_predict_test <- data.frame(xgb_CD8_TN_CN_predict_test)
# print(lm(xgb_CD8_TN_CN_predict_test$xgb_CD8_TN_CN_predict_test~CN_CD8_Tn_test$CN))
# print(cor(xgb_CD8_TN_CN_predict_test$xgb_CD8_TN_CN_predict_test,CN_CD8_Tn_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # R2 = 0.5532536
# 
# plot(xgb_CD8_TN_CN_predict_test$xgb_CD8_TN_CN_predict_test,CN_CD8_Tn_test$CN)
# 
# 
# 
# gc()
# rm(xgb_CD8_TN_CN)
# 
# 


plot(CN_CD8_Tn_with_RNA_train$TPM,CN_CD8_Tn_with_RNA_train$CN)
print(cor(CN_CD8_Tn_with_RNA_train$TPM,CN_CD8_Tn_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)








###__________________________________________________________________________________________###
###____________________________________________CD8_Tcm________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Tcm <- read.delim("~/Analysis/data/protein/CD8_Tcm_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tcm <- subset(CN_CD8_Tcm,CN_CD8_Tcm$CN>0)
dim(CN_CD8_Tcm)


TPM_CD8_Tcm <- read.delim("~/Analysis/data/RNA/CD8_Tcm_TPM_log10.csv",sep=";",dec=",")
TPM_CD8_Tcm$TPM <- 10^TPM_CD8_Tcm$TPM
TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD8_Tcm <- ddply(TPM_CD8_Tcm,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Tcm$TPM <- log10(TPM_CD8_Tcm$TPM)
#dim(TPM_CD8_Tcm[duplicated(TPM_CD8_Tcm$uniprotswissprot),])



CN_CD8_Tcm <- merge(CN_CD8_Tcm,TPM_CD8_Tcm,by.x="ID",by.y="uniprotswissprot")
CN_CD8_Tcm$ensembl_gene_id <- NULL
CN_CD8_Tcm <- merge(CN_CD8_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tcm[is.na(CN_CD8_Tcm)]= 0
print(dim(CN_CD8_Tcm)) # 7984 7126

CN_CD8_Tcm_dups <- CN_CD8_Tcm[duplicated(CN_CD8_Tcm$ID),]

rownames(CN_CD8_Tcm) <- CN_CD8_Tcm$ID
CN_CD8_Tcm$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tcm[is.na(CN_CD8_Tcm)]=0
CN_CD8_Tcm <- CN_CD8_Tcm[, -nearZeroVar(CN_CD8_Tcm, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD8_Tcm)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD8_Tcm), 0.8 * nrow(CN_CD8_Tcm))
CN_CD8_Tcm_train <- CN_CD8_Tcm[train_row,]
CN_CD8_Tcm_test <- CN_CD8_Tcm[-train_row,]

dim(CN_CD8_Tcm_train)
dim(CN_CD8_Tcm_test)

write.table(CN_CD8_Tcm_train,"~/Analysis/split_dataset/CN_CD8_Tcm_with_RNA_train",sep=";")
write.table(CN_CD8_Tcm_test,"~/Analysis/split_dataset/CN_CD8_Tcm_with_RNA_test",sep=";")

# CN_CD8_Tcm_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tcm_with_RNA_train",sep=";")
# CN_CD8_Tcm_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tcm_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
print("modeling CD8 Tcm CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8_Tcm_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                 max_depth = 6,
                                 colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                 ## The values below are default values in the sklearn-api.
                                 eta = 0.05,
                                 gamma=0,
                                 min_child_weight = 0.9,
                                 subsample = 1)

start_time <- Sys.time()
xgb_CD8_Tcm_CN <- train(CN~.,
                       data=CN_CD8_Tcm_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_CD8_Tcm_CN,
                       na.action = na.omit,
                       nthread=80,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD8_Tcm_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tcm_with_RNA_25_03_2022.RDS")
# xgb_CD8_Tcm_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tcm_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD8_Tcm_CN)

xgb_CD8_Tcm_CN_predict_test <- predict(xgb_CD8_Tcm_CN,CN_CD8_Tcm_test)
xgb_CD8_Tcm_CN_predict_test <- data.frame(xgb_CD8_Tcm_CN_predict_test)
print(lm(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test~CN_CD8_Tcm_test$CN))
print(cor(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test,CN_CD8_Tcm_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD8_Tcm_CN_predict_test$xgb_CD8_Tcm_CN_predict_test,CN_CD8_Tcm_test$CN)



gc()
rm(xgb_CD8_Tcm_CN)




# plot(CN_CD8_Tcm_with_RNA_train$TPM,CN_CD8_Tcm_with_RNA_train$CN)
# print(cor(CN_CD8_Tcm_with_RNA_train$TPM,CN_CD8_Tcm_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)











###__________________________________________________________________________________________###
###____________________________________________CD8_Tem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Tem <- read.delim("~/Analysis/data/protein/CD8_Tem_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tem <- subset(CN_CD8_Tem,CN_CD8_Tem$CN>0)
dim(CN_CD8_Tem)


TPM_CD8_Tem <- read.delim("~/Analysis/data/RNA/CD8_Tem_TPM_log10.csv",sep=";",dec=",")
TPM_CD8_Tem$TPM <- 10^TPM_CD8_Tem$TPM
TPM_CD8_Tem <- merge(TPM_CD8_Tem,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD8_Tem <- merge(TPM_CD8_Tem,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD8_Tem <- ddply(TPM_CD8_Tem,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Tem$TPM <- log10(TPM_CD8_Tem$TPM)
#dim(TPM_CD8_Tem[duplicated(TPM_CD8_Tem$uniprotswissprot),])



CN_CD8_Tem <- merge(CN_CD8_Tem,TPM_CD8_Tem,by.x="ID",by.y="uniprotswissprot")
CN_CD8_Tem$ensembl_gene_id <- NULL
CN_CD8_Tem <- merge(CN_CD8_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tem[is.na(CN_CD8_Tem)]= 0
print(dim(CN_CD8_Tem)) # 7984 7126

CN_CD8_Tem_dups <- CN_CD8_Tem[duplicated(CN_CD8_Tem$ID),]

rownames(CN_CD8_Tem) <- CN_CD8_Tem$ID
CN_CD8_Tem$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tem[is.na(CN_CD8_Tem)]=0
CN_CD8_Tem <- CN_CD8_Tem[, -nearZeroVar(CN_CD8_Tem, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD8_Tem)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD8_Tem), 0.8 * nrow(CN_CD8_Tem))
CN_CD8_Tem_train <- CN_CD8_Tem[train_row,]
CN_CD8_Tem_test <- CN_CD8_Tem[-train_row,]

dim(CN_CD8_Tem_train)
dim(CN_CD8_Tem_test)

write.table(CN_CD8_Tem_train,"~/Analysis/split_dataset/CN_CD8_Tem_with_RNA_train",sep=";")
write.table(CN_CD8_Tem_test,"~/Analysis/split_dataset/CN_CD8_Tem_with_RNA_test",sep=";")

# CN_CD8_Tem_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Tem_with_RNA_train",sep=";")
# CN_CD8_Tem_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Tem_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
print("modeling CD8 Tem CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8_Tem_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_CD8_Tem_CN <- train(CN~.,
                        data=CN_CD8_Tem_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD8_Tem_CN,
                        na.action = na.omit,
                        nthread=80,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD8_Tem_CN, "~/Analysis/models/protein/xgb_CN_CD8_Tem_with_RNA_25_03_2022.RDS")
# xgb_CD8_Tem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Tem_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD8_Tem_CN)

xgb_CD8_Tem_CN_predict_test <- predict(xgb_CD8_Tem_CN,CN_CD8_Tem_test)
xgb_CD8_Tem_CN_predict_test <- data.frame(xgb_CD8_Tem_CN_predict_test)
print(lm(xgb_CD8_Tem_CN_predict_test$xgb_CD8_Tem_CN_predict_test~CN_CD8_Tem_test$CN))
print(cor(xgb_CD8_Tem_CN_predict_test$xgb_CD8_Tem_CN_predict_test,CN_CD8_Tem_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD8_Tem_CN_predict_test$xgb_CD8_Tem_CN_predict_test,CN_CD8_Tem_test$CN)



gc()
rm(xgb_CD8_Tem_CN)


# plot(CN_CD8_Tem_with_RNA_train$TPM,CN_CD8_Tem_with_RNA_train$CN)
# print(cor(CN_CD8_Tem_with_RNA_train$TPM,CN_CD8_Tem_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)











###__________________________________________________________________________________________###
###____________________________________________CD8_Teff________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Teff <- read.delim("~/Analysis/data/protein/CD8_Teff_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Teff <- subset(CN_CD8_Teff,CN_CD8_Teff$CN>0)
dim(CN_CD8_Teff)


TPM_CD8_Teff <- read.delim("~/Analysis/data/RNA/CD8_Teff_TPM_log10.csv",sep=";",dec=",")
TPM_CD8_Teff$TPM <- 10^TPM_CD8_Teff$TPM
TPM_CD8_Teff <- merge(TPM_CD8_Teff,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD8_Teff <- merge(TPM_CD8_Teff,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD8_Teff <- ddply(TPM_CD8_Teff,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Teff$TPM <- log10(TPM_CD8_Teff$TPM)
#dim(TPM_CD8_Teff[duplicated(TPM_CD8_Teff$uniprotswissprot),])



CN_CD8_Teff <- merge(CN_CD8_Teff,TPM_CD8_Teff,by.x="ID",by.y="uniprotswissprot")
CN_CD8_Teff$ensembl_gene_id <- NULL
CN_CD8_Teff <- merge(CN_CD8_Teff,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Teff[is.na(CN_CD8_Teff)]= 0
print(dim(CN_CD8_Teff)) # 7984 7126

CN_CD8_Teff_dups <- CN_CD8_Teff[duplicated(CN_CD8_Teff$ID),]

rownames(CN_CD8_Teff) <- CN_CD8_Teff$ID
CN_CD8_Teff$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Teff[is.na(CN_CD8_Teff)]=0
CN_CD8_Teff <- CN_CD8_Teff[, -nearZeroVar(CN_CD8_Teff, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_CD8_Teff)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_CD8_Teff), 0.8 * nrow(CN_CD8_Teff))
CN_CD8_Teff_train <- CN_CD8_Teff[train_row,]
CN_CD8_Teff_test <- CN_CD8_Teff[-train_row,]

dim(CN_CD8_Teff_train)
dim(CN_CD8_Teff_test)

write.table(CN_CD8_Teff_train,"~/Analysis/split_dataset/CN_CD8_Teff_with_RNA_train",sep=";")
write.table(CN_CD8_Teff_test,"~/Analysis/split_dataset/CN_CD8_Teff_with_RNA_test",sep=";")

# CN_CD8_Teff_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD8_Teff_with_RNA_train",sep=";")
# CN_CD8_Teff_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD8_Teff_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
print("modeling CD8 Teff CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_CD8_Teff_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_CD8_Teff_CN <- train(CN~.,
                        data=CN_CD8_Teff_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_CD8_Teff_CN,
                        na.action = na.omit,
                        nthread=80,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD8_Teff_CN, "~/Analysis/models/protein/xgb_CN_CD8_Teff_with_RNA_25_03_2022.RDS")
# xgb_CD8_Teff_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD8_Teff_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD8_Teff_CN)

xgb_CD8_Teff_CN_predict_test <- predict(xgb_CD8_Teff_CN,CN_CD8_Teff_test)
xgb_CD8_Teff_CN_predict_test <- data.frame(xgb_CD8_Teff_CN_predict_test)
print(lm(xgb_CD8_Teff_CN_predict_test$xgb_CD8_Teff_CN_predict_test~CN_CD8_Teff_test$CN))
print(cor(xgb_CD8_Teff_CN_predict_test$xgb_CD8_Teff_CN_predict_test,CN_CD8_Teff_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD8_Teff_CN_predict_test$xgb_CD8_Teff_CN_predict_test,CN_CD8_Teff_test$CN)



gc()
rm(xgb_CD8_Teff_CN)





# plot(CN_CD8_Teff_with_RNA_train$TPM,CN_CD8_Teff_with_RNA_train$CN)
# print(cor(CN_CD8_Teff_with_RNA_train$TPM,CN_CD8_Teff_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)











###__________________________________________________________________________________________###
###____________________________________________CD4_Tn________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tn <- read.delim("~/Analysis/data/protein/CD4_Tn_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tn <- subset(CN_CD4_Tn,CN_CD4_Tn$CN>0)
dim(CN_CD4_Tn)


TPM_CD4_Tn <- read.delim("~/Analysis/data/RNA/CD4_Tn_TPM_log10.csv",sep=";",dec=",")
TPM_CD4_Tn$TPM <- 10^TPM_CD4_Tn$TPM
TPM_CD4_Tn <- merge(TPM_CD4_Tn,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD4_Tn <- merge(TPM_CD4_Tn,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD4_Tn <- ddply(TPM_CD4_Tn,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD4_Tn$TPM <- log10(TPM_CD4_Tn$TPM)
#dim(TPM_CD4_Tn[duplicated(TPM_CD4_Tn$uniprotswissprot),])



CN_CD4_Tn <- merge(CN_CD4_Tn,TPM_CD4_Tn,by.x="ID",by.y="uniprotswissprot")
CN_CD4_Tn$ensembl_gene_id <- NULL
CN_CD4_Tn <- merge(CN_CD4_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tn[is.na(CN_CD4_Tn)]= 0
print(dim(CN_CD4_Tn)) # 7984 7126

CN_CD4_Tn_dups <- CN_CD4_Tn[duplicated(CN_CD4_Tn$ID),]

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

write.table(CN_CD4_Tn_train,"~/Analysis/split_dataset/CN_CD4_Tn_with_RNA_train",sep=";")
write.table(CN_CD4_Tn_test,"~/Analysis/split_dataset/CN_CD4_Tn_with_RNA_test",sep=";")

# CN_CD4_Tn_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_with_RNA_train",sep=";")
# CN_CD4_Tn_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
print("modeling CD4 Tn CN")
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
                        nthread=80,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD4_Tn_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tn_with_RNA_25_03_2022.RDS")
# xgb_CD4_Tn_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tn_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD4_Tn_CN)

xgb_CD4_Tn_CN_predict_test <- predict(xgb_CD4_Tn_CN,CN_CD4_Tn_test)
xgb_CD4_Tn_CN_predict_test <- data.frame(xgb_CD4_Tn_CN_predict_test)
print(lm(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test~CN_CD4_Tn_test$CN))
print(cor(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN)



gc()
rm(xgb_CD4_Tn_CN)



# plot(CN_CD4_Tn_with_RNA_train$TPM,CN_CD4_Tn_with_RNA_train$CN)
# print(cor(CN_CD4_Tn_with_RNA_train$TPM,CN_CD4_Tn_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)











###__________________________________________________________________________________________###
###____________________________________________CD4_Tcm________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tcm <- read.delim("~/Analysis/data/protein/CD4_Tcm_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tcm <- subset(CN_CD4_Tcm,CN_CD4_Tcm$CN>0)
dim(CN_CD4_Tcm)


TPM_CD4_Tcm <- read.delim("~/Analysis/data/RNA/CD4_Tcm_TPM_log10.csv",sep=";",dec=",")
TPM_CD4_Tcm$TPM <- 10^TPM_CD4_Tcm$TPM
TPM_CD4_Tcm <- merge(TPM_CD4_Tcm,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD4_Tcm <- merge(TPM_CD4_Tcm,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD4_Tcm <- ddply(TPM_CD4_Tcm,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD4_Tcm$TPM <- log10(TPM_CD4_Tcm$TPM)
#dim(TPM_CD4_Tcm[duplicated(TPM_CD4_Tcm$uniprotswissprot),])



CN_CD4_Tcm <- merge(CN_CD4_Tcm,TPM_CD4_Tcm,by.x="ID",by.y="uniprotswissprot")
CN_CD4_Tcm$ensembl_gene_id <- NULL
CN_CD4_Tcm <- merge(CN_CD4_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tcm[is.na(CN_CD4_Tcm)]= 0
print(dim(CN_CD4_Tcm)) # 7984 7126

CN_CD4_Tcm_dups <- CN_CD4_Tcm[duplicated(CN_CD4_Tcm$ID),]

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

write.table(CN_CD4_Tcm_train,"~/Analysis/split_dataset/CN_CD4_Tcm_with_RNA_train",sep=";")
write.table(CN_CD4_Tcm_test,"~/Analysis/split_dataset/CN_CD4_Tcm_with_RNA_test",sep=";")

# CN_CD4_Tcm_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tcm_with_RNA_train",sep=";")
# CN_CD4_Tcm_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tcm_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
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
                       nthread=80,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD4_Tcm_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tcm_with_RNA_25_03_2022.RDS")
# xgb_CD4_Tcm_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tcm_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD4_Tcm_CN)

xgb_CD4_Tcm_CN_predict_test <- predict(xgb_CD4_Tcm_CN,CN_CD4_Tcm_test)
xgb_CD4_Tcm_CN_predict_test <- data.frame(xgb_CD4_Tcm_CN_predict_test)
print(lm(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test~CN_CD4_Tcm_test$CN))
print(cor(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test,CN_CD4_Tcm_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD4_Tcm_CN_predict_test$xgb_CD4_Tcm_CN_predict_test,CN_CD4_Tcm_test$CN)



gc()
rm(xgb_CD4_Tcm_CN)




# plot(CN_CD4_Tcm_with_RNA_train$TPM,CN_CD4_Tcm_with_RNA_train$CN)
# print(cor(CN_CD4_Tcm_with_RNA_train$TPM,CN_CD4_Tcm_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)









###__________________________________________________________________________________________###
###____________________________________________CD4_Tem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tem <- read.delim("~/Analysis/data/protein/CD4_Tem_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tem <- subset(CN_CD4_Tem,CN_CD4_Tem$CN>0)
dim(CN_CD4_Tem)


TPM_CD4_Tem <- read.delim("~/Analysis/data/RNA/CD4_Tem_TPM_log10.csv",sep=";",dec=",")
TPM_CD4_Tem$TPM <- 10^TPM_CD4_Tem$TPM
TPM_CD4_Tem <- merge(TPM_CD4_Tem,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_CD4_Tem <- merge(TPM_CD4_Tem,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_CD4_Tem <- ddply(TPM_CD4_Tem,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD4_Tem$TPM <- log10(TPM_CD4_Tem$TPM)
#dim(TPM_CD4_Tem[duplicated(TPM_CD4_Tem$uniprotswissprot),])



CN_CD4_Tem <- merge(CN_CD4_Tem,TPM_CD4_Tem,by.x="ID",by.y="uniprotswissprot")
CN_CD4_Tem$ensembl_gene_id <- NULL
CN_CD4_Tem <- merge(CN_CD4_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tem[is.na(CN_CD4_Tem)]= 0
print(dim(CN_CD4_Tem)) # 7984 7126

CN_CD4_Tem_dups <- CN_CD4_Tem[duplicated(CN_CD4_Tem$ID),]

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

write.table(CN_CD4_Tem_train,"~/Analysis/split_dataset/CN_CD4_Tem_with_RNA_train",sep=";")
write.table(CN_CD4_Tem_test,"~/Analysis/split_dataset/CN_CD4_Tem_with_RNA_test",sep=";")

# CN_CD4_Tem_with_RNA_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_with_RNA_train",sep=";")
# CN_CD4_Tem_with_RNA_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_with_RNA_test",sep=";")
# 
# 
# 
# ## Training the xgb model ##
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
                       nthread=80,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_CD4_Tem_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tem_with_RNA_25_03_2022.RDS")
# xgb_CD4_Tem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tem_with_RNA_25_03_2022.RDS")
# 
# print(xgb_CD4_Tem_CN)

xgb_CD4_Tem_CN_predict_test <- predict(xgb_CD4_Tem_CN,CN_CD4_Tem_test)
xgb_CD4_Tem_CN_predict_test <- data.frame(xgb_CD4_Tem_CN_predict_test)
print(lm(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test~CN_CD4_Tem_test$CN))
print(cor(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test,CN_CD4_Tem_test$CN, method = "pearson", use = "complete.obs")^2)

plot(xgb_CD4_Tem_CN_predict_test$xgb_CD4_Tem_CN_predict_test,CN_CD4_Tem_test$CN)



gc()
rm(xgb_CD4_Tem_CN)




plot(CN_CD4_Tem_with_RNA_train$TPM,CN_CD4_Tem_with_RNA_train$CN)
print(cor(CN_CD4_Tem_with_RNA_train$TPM,CN_CD4_Tem_with_RNA_train$CN, method = "pearson", use = "complete.obs")^2)






