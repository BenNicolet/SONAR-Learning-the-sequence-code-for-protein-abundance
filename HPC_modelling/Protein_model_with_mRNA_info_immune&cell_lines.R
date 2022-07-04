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
###____________________________________________Bnaive________________________________________###
###__________________________________________________________________________________________###


# ## prep data ##
# CN_Bnaive <- read.delim("~/Analysis/data/protein/B_naive_CN_log10.csv",sep=";",dec = ",")
# CN_Bnaive <- subset(CN_Bnaive,CN_Bnaive$CN>0)
# dim(CN_Bnaive)
# 
# 
# TPM_Bnaive <- read.delim("~/Analysis/data/immune_cells_RNA/B_naive_TPM_log10.csv",sep=";",dec=",")
# TPM_Bnaive$TPM <- 10^TPM_Bnaive$TPM
# TPM_Bnaive <- merge(TPM_Bnaive,tx2gene, by.x="ID", by.y="ensembl_transcript_id")
# 
# 
# registerDoMC(10)
# 
# TPM_Bnaive <- merge(TPM_Bnaive,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
# TPM_Bnaive <- ddply(TPM_Bnaive,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Bnaive$TPM <- log10(TPM_Bnaive$TPM)
# #dim(TPM_Bnaive[duplicated(TPM_Bnaive$uniprotswissprot),])
# 
# 
# 
# CN_Bnaive <- merge(CN_Bnaive,TPM_Bnaive,by.x="ID",by.y="uniprotswissprot")
# CN_Bnaive$ensembl_gene_id <- NULL
# CN_Bnaive <- merge(CN_Bnaive,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
# CN_Bnaive[is.na(CN_Bnaive)]= 0
# print(dim(CN_Bnaive)) # 7984 7126
# 
# CN_Bnaive_dups <- CN_Bnaive[duplicated(CN_Bnaive$ID),]
# 
# rownames(CN_Bnaive) <- CN_Bnaive$ID
# CN_Bnaive$ID <- NULL
# 
# 
# ## Feature selection ##
# registerDoMC(6)
# 
# CN_Bnaive[is.na(CN_Bnaive)]=0
# CN_Bnaive <- CN_Bnaive[, -nearZeroVar(CN_Bnaive, allowParallel = T, uniqueCut = 0.1)]
# print(dim(CN_Bnaive)) # 7984 2759
# 
# gc()
# 
# 
# 
# 
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(CN_Bnaive), 0.8 * nrow(CN_Bnaive))
# CN_Bnaive_train <- CN_Bnaive[train_row,]
# CN_Bnaive_test <- CN_Bnaive[-train_row,]
# 
# dim(CN_Bnaive_train)
# dim(CN_Bnaive_test)
# 
# write.table(CN_Bnaive_train,"~/Analysis/split_dataset/CN_Bnaive_with_RNA_train",sep=";")
# write.table(CN_Bnaive_test,"~/Analysis/split_dataset/CN_Bnaive_with_RNA_test",sep=";")
# 
# # CN_Bnaive_train <- read.delim("~/Analysis/split_dataset/CN_Bnaive_with_RNA_train",sep=";")
# # CN_Bnaive_test <- read.delim("~/Analysis/split_dataset/CN_Bnaive_with_RNA_test",sep=";")
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
# xgbGrid_Bnaive_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                               max_depth = 6,
#                               colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                               ## The values below are default values in the sklearn-api.
#                               eta = 0.05,
#                               gamma=0,
#                               min_child_weight = 0.9,
#                               subsample = 1)
# 
# start_time <- Sys.time()
# xgb_Bnaive_CN <- train(CN~.,
#                         data=CN_Bnaive_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_Bnaive_CN,
#                         na.action = na.omit,
#                         nthread=90,
#                         verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# 
# saveRDS(xgb_Bnaive_CN, "~/Analysis/models/protein/xgb_CN_Bnaive_with_RNA_02_06_2022.RDS")
# # xgb_Bnaive_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Bnaive_with_RNA_02_06_2022.RDS")
# #
# # print(xgb_Bnaive_CN)
# 
# xgb_Bnaive_CN_predict_test <- predict(xgb_Bnaive_CN,CN_Bnaive_test)
# xgb_Bnaive_CN_predict_test <- data.frame(xgb_Bnaive_CN_predict_test)
# # print(lm(xgb_Bnaive_CN_predict_test$xgb_Bnaive_CN_predict_test~CN_Bnaive_test$CN))
# print(cor(xgb_Bnaive_CN_predict_test$xgb_Bnaive_CN_predict_test,CN_Bnaive_test$CN, method = "pearson", use = "complete.obs")^2)
# 
# # plot(CN_Bnaive_train$TPM,CN_Bnaive_train$CN)
# print(cor(CN_Bnaive_train$TPM,CN_Bnaive_train$CN, method = "pearson", use = "complete.obs")^2)
# 
plot(xgb_Bnaive_CN_predict_test$xgb_Bnaive_CN_predict_test,CN_Bnaive_test$CN,asp=1)
# 
# gc()
# rm(xgb_Bnaive_CN)
















###__________________________________________________________________________________________###
###____________________________________________B_mem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_B_mem <- read.delim("~/Analysis/data/protein/B_mem_CN_log10.csv",sep=";",dec = ",")
CN_B_mem <- subset(CN_B_mem,CN_B_mem$CN>0)
dim(CN_B_mem)


TPM_B_NSmem <- read.delim("~/Analysis/data/immune_cells_RNA/B_NS_mem_TPM_log10.csv",sep=";",dec=",")
TPM_B_Smem <- read.delim("~/Analysis/data/immune_cells_RNA/B_S_mem_TPM_log10.csv",sep=";",dec=",")
TPM_B_mem <- merge(TPM_B_NSmem,TPM_B_Smem,by="ID")

TPM_B_mem$TPM <- (10^TPM_B_mem$TPM.x + 10^TPM_B_mem$TPM.y)/2
TPM_B_mem$TPM.x <- NULL
TPM_B_mem$TPM.y <- NULL

TPM_B_mem <- merge(TPM_B_mem,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_B_mem <- merge(TPM_B_mem,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_B_mem <- ddply(TPM_B_mem,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_B_mem$TPM <- log10(TPM_B_mem$TPM)
#dim(TPM_B_mem[duplicated(TPM_B_mem$uniprotswissprot),])



CN_B_mem <- merge(CN_B_mem,TPM_B_mem,by.x="ID",by.y="uniprotswissprot")
CN_B_mem$ensembl_gene_id <- NULL
CN_B_mem <- merge(CN_B_mem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_B_mem[is.na(CN_B_mem)]= 0
print(dim(CN_B_mem)) # 7984 7126

CN_B_mem_dups <- CN_B_mem[duplicated(CN_B_mem$ID),]

rownames(CN_B_mem) <- CN_B_mem$ID
CN_B_mem$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_B_mem[is.na(CN_B_mem)]=0
CN_B_mem <- CN_B_mem[, -nearZeroVar(CN_B_mem, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_B_mem)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_B_mem), 0.8 * nrow(CN_B_mem))
CN_B_mem_train <- CN_B_mem[train_row,]
CN_B_mem_test <- CN_B_mem[-train_row,]

dim(CN_B_mem_train)
dim(CN_B_mem_test)

write.table(CN_B_mem_train,"~/Analysis/split_dataset/CN_B_mem_with_RNA_train",sep=";")
write.table(CN_B_mem_test,"~/Analysis/split_dataset/CN_B_mem_with_RNA_test",sep=";")

# CN_B_mem_train <- read.delim("~/Analysis/split_dataset/CN_B_mem_with_RNA_train",sep=";")
# CN_B_mem_test <- read.delim("~/Analysis/split_dataset/CN_B_mem_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_B_mem_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                 max_depth = 6,
                                 colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                 ## The values below are default values in the sklearn-api.
                                 eta = 0.05,
                                 gamma=0,
                                 min_child_weight = 0.9,
                                 subsample = 1)

start_time <- Sys.time()
xgb_B_mem_CN <- train(CN~.,
                       data=CN_B_mem_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_B_mem_CN,
                       na.action = na.omit,
                       nthread=90,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_B_mem_CN, "~/Analysis/models/protein/xgb_CN_B_mem_with_RNA_02_06_2022.RDS")
# xgb_B_mem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_B_mem_with_RNA_02_06_2022.RDS")
#
# print(xgb_B_mem_CN)

xgb_B_mem_CN_predict_test <- predict(xgb_B_mem_CN,CN_B_mem_test)
xgb_B_mem_CN_predict_test <- data.frame(xgb_B_mem_CN_predict_test)
# print(lm(xgb_B_mem_CN_predict_test$xgb_B_mem_CN_predict_test~CN_B_mem_test$CN))
print(cor(xgb_B_mem_CN_predict_test$xgb_B_mem_CN_predict_test,CN_B_mem_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_B_mem_train$TPM,CN_B_mem_train$CN)
# hist(CN_B_mem_train$CN)
print(cor(CN_B_mem_train$TPM,CN_B_mem_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_B_mem_CN_predict_test$xgb_B_mem_CN_predict_test,CN_B_mem_test$CN)

gc()
rm(xgb_B_mem_CN)












###__________________________________________________________________________________________###
###____________________________________________B_plasma________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_B_plasma <- read.delim("~/Analysis/data/protein/B_plasma_CN_log10.csv",sep=";",dec = ",")
CN_B_plasma <- subset(CN_B_plasma,CN_B_plasma$CN>0)
dim(CN_B_plasma)


TPM_B_plasma <- read.delim("~/Analysis/data/immune_cells_RNA/B_plasma_TPM_log10.csv",sep=";",dec=",")
TPM_B_plasma$TPM <- 10^TPM_B_plasma$TPM
TPM_B_plasma <- merge(TPM_B_plasma,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_B_plasma <- merge(TPM_B_plasma,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_B_plasma <- ddply(TPM_B_plasma,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_B_plasma$TPM <- log10(TPM_B_plasma$TPM)
#dim(TPM_B_plasma[duplicated(TPM_B_plasma$uniprotswissprot),])



CN_B_plasma <- merge(CN_B_plasma,TPM_B_plasma,by.x="ID",by.y="uniprotswissprot")
CN_B_plasma$ensembl_gene_id <- NULL
CN_B_plasma <- merge(CN_B_plasma,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_B_plasma[is.na(CN_B_plasma)]= 0
print(dim(CN_B_plasma)) # 7984 7126

CN_B_plasma_dups <- CN_B_plasma[duplicated(CN_B_plasma$ID),]

rownames(CN_B_plasma) <- CN_B_plasma$ID
CN_B_plasma$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_B_plasma[is.na(CN_B_plasma)]=0
CN_B_plasma <- CN_B_plasma[, -nearZeroVar(CN_B_plasma, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_B_plasma)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_B_plasma), 0.8 * nrow(CN_B_plasma))
CN_B_plasma_train <- CN_B_plasma[train_row,]
CN_B_plasma_test <- CN_B_plasma[-train_row,]

dim(CN_B_plasma_train)
dim(CN_B_plasma_test)

write.table(CN_B_plasma_train,"~/Analysis/split_dataset/CN_B_plasma_with_RNA_train",sep=";")
write.table(CN_B_plasma_test,"~/Analysis/split_dataset/CN_B_plasma_with_RNA_test",sep=";")

# CN_B_plasma_train <- read.delim("~/Analysis/split_dataset/CN_B_plasma_with_RNA_train",sep=";")
# CN_B_plasma_test <- read.delim("~/Analysis/split_dataset/CN_B_plasma_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_B_plasma_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                 max_depth = 6,
                                 colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                 ## The values below are default values in the sklearn-api.
                                 eta = 0.05,
                                 gamma=0,
                                 min_child_weight = 0.9,
                                 subsample = 1)

start_time <- Sys.time()
xgb_B_plasma_CN <- train(CN~.,
                       data=CN_B_plasma_train,
                       method="xgbTree",
                       trControl=control,
                       #metric="Rsquared",
                       tuneGrid= xgbGrid_B_plasma_CN,
                       na.action = na.omit,
                       nthread=90,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_B_plasma_CN, "~/Analysis/models/protein/xgb_CN_B_plasma_with_RNA_02_06_2022.RDS")
# xgb_B_plasma_CN <- readRDS("~/Analysis/models/protein/xgb_CN_B_plasma_with_RNA_02_06_2022.RDS")
#
# print(xgb_B_plasma_CN)

xgb_B_plasma_CN_predict_test <- predict(xgb_B_plasma_CN,CN_B_plasma_test)
xgb_B_plasma_CN_predict_test <- data.frame(xgb_B_plasma_CN_predict_test)
# print(lm(xgb_B_plasma_CN_predict_test$xgb_B_plasma_CN_predict_test~CN_B_plasma_test$CN))
print(cor(xgb_B_plasma_CN_predict_test$xgb_B_plasma_CN_predict_test,CN_B_plasma_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_B_plasma_with_RNA_train$TPM,CN_B_plasma_with_RNA_train$CN)
print(cor(CN_B_plasma_train$TPM,CN_B_plasma_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_B_plasma_CN_predict_test$xgb_B_plasma_CN_predict_test,CN_B_plasma_test$CN)

gc()
rm(xgb_B_plasma_CN)

















###__________________________________________________________________________________________###
###____________________________________________Basophils________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Basophils <- read.delim("~/Analysis/data/protein/Basophils_CN_log10.csv",sep=";",dec = ",")
CN_Basophils <- subset(CN_Basophils,CN_Basophils$CN>0)
dim(CN_Basophils)


TPM_Basophils <- read.delim("~/Analysis/data/immune_cells_RNA/Basophils_TPM_log10.csv",sep=";",dec=",")
TPM_Basophils$TPM <- 10^TPM_Basophils$TPM
TPM_Basophils <- merge(TPM_Basophils,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Basophils <- merge(TPM_Basophils,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Basophils <- ddply(TPM_Basophils,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Basophils$TPM <- log10(TPM_Basophils$TPM)
#dim(TPM_Basophils[duplicated(TPM_Basophils$uniprotswissprot),])



CN_Basophils <- merge(CN_Basophils,TPM_Basophils,by.x="ID",by.y="uniprotswissprot")
CN_Basophils$ensembl_gene_id <- NULL
CN_Basophils <- merge(CN_Basophils,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Basophils[is.na(CN_Basophils)]= 0
print(dim(CN_Basophils)) # 7984 7126

CN_Basophils_dups <- CN_Basophils[duplicated(CN_Basophils$ID),]

rownames(CN_Basophils) <- CN_Basophils$ID
CN_Basophils$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Basophils[is.na(CN_Basophils)]=0
CN_Basophils <- CN_Basophils[, -nearZeroVar(CN_Basophils, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Basophils)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Basophils), 0.8 * nrow(CN_Basophils))
CN_Basophils_train <- CN_Basophils[train_row,]
CN_Basophils_test <- CN_Basophils[-train_row,]

dim(CN_Basophils_train)
dim(CN_Basophils_test)

write.table(CN_Basophils_train,"~/Analysis/split_dataset/CN_Basophils_with_RNA_train",sep=";")
write.table(CN_Basophils_test,"~/Analysis/split_dataset/CN_Basophils_with_RNA_test",sep=";")

# CN_Basophils_train <- read.delim("~/Analysis/split_dataset/CN_Basophils_with_RNA_train",sep=";")
# CN_Basophils_test <- read.delim("~/Analysis/split_dataset/CN_Basophils_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Basophils_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_Basophils_CN <- train(CN~.,
                         data=CN_Basophils_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_Basophils_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Basophils_CN, "~/Analysis/models/protein/xgb_CN_Basophils_with_RNA_02_06_2022.RDS")
# xgb_Basophils_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Basophils_with_RNA_02_06_2022.RDS")
#
# print(xgb_Basophils_CN)

xgb_Basophils_CN_predict_test <- predict(xgb_Basophils_CN,CN_Basophils_test)
xgb_Basophils_CN_predict_test <- data.frame(xgb_Basophils_CN_predict_test)
# print(lm(xgb_Basophils_CN_predict_test$xgb_Basophils_CN_predict_test~CN_Basophils_test$CN))
print(cor(xgb_Basophils_CN_predict_test$xgb_Basophils_CN_predict_test,CN_Basophils_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Basophils_with_RNA_train$TPM,CN_Basophils_with_RNA_train$CN)
print(cor(CN_Basophils_train$TPM,CN_Basophils_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Basophils_CN_predict_test$xgb_Basophils_CN_predict_test,CN_Basophils_test$CN)

gc()
rm(xgb_Basophils_CN)


















###__________________________________________________________________________________________###
###____________________________________________Neutrophils________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Neutrophils <- read.delim("~/Analysis/data/protein/Neutrophils_CN_log10.csv",sep=";",dec = ",")
CN_Neutrophils <- subset(CN_Neutrophils,CN_Neutrophils$CN>0)
dim(CN_Neutrophils)


TPM_Neutrophils <- read.delim("~/Analysis/data/immune_cells_RNA/Neutrophils_TPM_log10.csv",sep=";",dec=",")
TPM_Neutrophils$TPM <- 10^TPM_Neutrophils$TPM
TPM_Neutrophils <- merge(TPM_Neutrophils,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Neutrophils <- merge(TPM_Neutrophils,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Neutrophils <- ddply(TPM_Neutrophils,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Neutrophils$TPM <- log10(TPM_Neutrophils$TPM)
#dim(TPM_Neutrophils[duplicated(TPM_Neutrophils$uniprotswissprot),])



CN_Neutrophils <- merge(CN_Neutrophils,TPM_Neutrophils,by.x="ID",by.y="uniprotswissprot")
CN_Neutrophils$ensembl_gene_id <- NULL
CN_Neutrophils <- merge(CN_Neutrophils,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Neutrophils[is.na(CN_Neutrophils)]= 0
print(dim(CN_Neutrophils)) # 7984 7126

CN_Neutrophils_dups <- CN_Neutrophils[duplicated(CN_Neutrophils$ID),]

rownames(CN_Neutrophils) <- CN_Neutrophils$ID
CN_Neutrophils$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Neutrophils[is.na(CN_Neutrophils)]=0
CN_Neutrophils <- CN_Neutrophils[, -nearZeroVar(CN_Neutrophils, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Neutrophils)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Neutrophils), 0.8 * nrow(CN_Neutrophils))
CN_Neutrophils_train <- CN_Neutrophils[train_row,]
CN_Neutrophils_test <- CN_Neutrophils[-train_row,]

dim(CN_Neutrophils_train)
dim(CN_Neutrophils_test)

write.table(CN_Neutrophils_train,"~/Analysis/split_dataset/CN_Neutrophils_with_RNA_train",sep=";")
write.table(CN_Neutrophils_test,"~/Analysis/split_dataset/CN_Neutrophils_with_RNA_test",sep=";")

# CN_Neutrophils_train <- read.delim("~/Analysis/split_dataset/CN_Neutrophils_with_RNA_train",sep=";")
# CN_Neutrophils_test <- read.delim("~/Analysis/split_dataset/CN_Neutrophils_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Neutrophils_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_Neutrophils_CN <- train(CN~.,
                         data=CN_Neutrophils_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_Neutrophils_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Neutrophils_CN, "~/Analysis/models/protein/xgb_CN_Neutrophils_with_RNA_02_06_2022.RDS")
# xgb_Neutrophils_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Neutrophils_with_RNA_02_06_2022.RDS")
#
# print(xgb_Neutrophils_CN)

xgb_Neutrophils_CN_predict_test <- predict(xgb_Neutrophils_CN,CN_Neutrophils_test)
xgb_Neutrophils_CN_predict_test <- data.frame(xgb_Neutrophils_CN_predict_test)
# print(lm(xgb_Neutrophils_CN_predict_test$xgb_Neutrophils_CN_predict_test~CN_Neutrophils_test$CN))
print(cor(xgb_Neutrophils_CN_predict_test$xgb_Neutrophils_CN_predict_test,CN_Neutrophils_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Neutrophils_with_RNA_train$TPM,CN_Neutrophils_with_RNA_train$CN)
print(cor(CN_Neutrophils_train$TPM,CN_Neutrophils_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Neutrophils_CN_predict_test$xgb_Neutrophils_CN_predict_test,CN_Neutrophils_test$CN)

gc()
rm(xgb_Neutrophils_CN)
















###__________________________________________________________________________________________###
###____________________________________________mDC________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_mDC <- read.delim("~/Analysis/data/protein/mDC_CN_log10.csv",sep=";",dec = ",")
CN_mDC <- subset(CN_mDC,CN_mDC$CN>0)
dim(CN_mDC)


TPM_mDC <- read.delim("~/Analysis/data/immune_cells_RNA/mDC_TPM_log10.csv",sep=";",dec=",")
TPM_mDC$TPM <- 10^TPM_mDC$TPM
TPM_mDC <- merge(TPM_mDC,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_mDC <- merge(TPM_mDC,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_mDC <- ddply(TPM_mDC,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_mDC$TPM <- log10(TPM_mDC$TPM)
#dim(TPM_mDC[duplicated(TPM_mDC$uniprotswissprot),])



CN_mDC <- merge(CN_mDC,TPM_mDC,by.x="ID",by.y="uniprotswissprot")
CN_mDC$ensembl_gene_id <- NULL
CN_mDC <- merge(CN_mDC,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_mDC[is.na(CN_mDC)]= 0
print(dim(CN_mDC)) # 7984 7126

CN_mDC_dups <- CN_mDC[duplicated(CN_mDC$ID),]

rownames(CN_mDC) <- CN_mDC$ID
CN_mDC$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_mDC[is.na(CN_mDC)]=0
CN_mDC <- CN_mDC[, -nearZeroVar(CN_mDC, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_mDC)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_mDC), 0.8 * nrow(CN_mDC))
CN_mDC_train <- CN_mDC[train_row,]
CN_mDC_test <- CN_mDC[-train_row,]

dim(CN_mDC_train)
dim(CN_mDC_test)

write.table(CN_mDC_train,"~/Analysis/split_dataset/CN_mDC_with_RNA_train",sep=";")
write.table(CN_mDC_test,"~/Analysis/split_dataset/CN_mDC_with_RNA_test",sep=";")

# CN_mDC_train <- read.delim("~/Analysis/split_dataset/CN_mDC_with_RNA_train",sep=";")
# CN_mDC_test <- read.delim("~/Analysis/split_dataset/CN_mDC_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_mDC_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_mDC_CN <- train(CN~.,
                         data=CN_mDC_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_mDC_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_mDC_CN, "~/Analysis/models/protein/xgb_CN_mDC_with_RNA_02_06_2022.RDS")
# xgb_mDC_CN <- readRDS("~/Analysis/models/protein/xgb_CN_mDC_with_RNA_02_06_2022.RDS")
#
# print(xgb_mDC_CN)

xgb_mDC_CN_predict_test <- predict(xgb_mDC_CN,CN_mDC_test)
xgb_mDC_CN_predict_test <- data.frame(xgb_mDC_CN_predict_test)
# print(lm(xgb_mDC_CN_predict_test$xgb_mDC_CN_predict_test~CN_mDC_test$CN))
print(cor(xgb_mDC_CN_predict_test$xgb_mDC_CN_predict_test,CN_mDC_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_mDC_with_RNA_train$TPM,CN_mDC_with_RNA_train$CN)
print(cor(CN_mDC_train$TPM,CN_mDC_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_mDC_CN_predict_test$xgb_mDC_CN_predict_test,CN_mDC_test$CN)

gc()
rm(xgb_mDC_CN)
















###__________________________________________________________________________________________###
###____________________________________________pDC________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_pDC <- read.delim("~/Analysis/data/protein/pDC_CN_log10.csv",sep=";",dec = ",")
CN_pDC <- subset(CN_pDC,CN_pDC$CN>0)
dim(CN_pDC)


TPM_pDC <- read.delim("~/Analysis/data/immune_cells_RNA/pDC_TPM_log10.csv",sep=";",dec=",")
TPM_pDC$TPM <- 10^TPM_pDC$TPM
TPM_pDC <- merge(TPM_pDC,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_pDC <- merge(TPM_pDC,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_pDC <- ddply(TPM_pDC,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_pDC$TPM <- log10(TPM_pDC$TPM)
#dim(TPM_pDC[duplicated(TPM_pDC$uniprotswissprot),])



CN_pDC <- merge(CN_pDC,TPM_pDC,by.x="ID",by.y="uniprotswissprot")
CN_pDC$ensembl_gene_id <- NULL
CN_pDC <- merge(CN_pDC,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_pDC[is.na(CN_pDC)]= 0
print(dim(CN_pDC)) # 7984 7126

CN_pDC_dups <- CN_pDC[duplicated(CN_pDC$ID),]

rownames(CN_pDC) <- CN_pDC$ID
CN_pDC$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_pDC[is.na(CN_pDC)]=0
CN_pDC <- CN_pDC[, -nearZeroVar(CN_pDC, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_pDC)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_pDC), 0.8 * nrow(CN_pDC))
CN_pDC_train <- CN_pDC[train_row,]
CN_pDC_test <- CN_pDC[-train_row,]

dim(CN_pDC_train)
dim(CN_pDC_test)

write.table(CN_pDC_train,"~/Analysis/split_dataset/CN_pDC_with_RNA_train",sep=";")
write.table(CN_pDC_test,"~/Analysis/split_dataset/CN_pDC_with_RNA_test",sep=";")

# CN_pDC_train <- read.delim("~/Analysis/split_dataset/CN_pDC_with_RNA_train",sep=";")
# CN_pDC_test <- read.delim("~/Analysis/split_dataset/CN_pDC_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_pDC_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_pDC_CN <- train(CN~.,
                         data=CN_pDC_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_pDC_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_pDC_CN, "~/Analysis/models/protein/xgb_CN_pDC_with_RNA_02_06_2022.RDS")
# xgb_pDC_CN <- readRDS("~/Analysis/models/protein/xgb_CN_pDC_with_RNA_02_06_2022.RDS")
#
# print(xgb_pDC_CN)

xgb_pDC_CN_predict_test <- predict(xgb_pDC_CN,CN_pDC_test)
xgb_pDC_CN_predict_test <- data.frame(xgb_pDC_CN_predict_test)
# print(lm(xgb_pDC_CN_predict_test$xgb_pDC_CN_predict_test~CN_pDC_test$CN))
print(cor(xgb_pDC_CN_predict_test$xgb_pDC_CN_predict_test,CN_pDC_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_pDC_with_RNA_train$TPM,CN_pDC_with_RNA_train$CN)
print(cor(CN_pDC_train$TPM,CN_pDC_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_pDC_CN_predict_test$xgb_pDC_CN_predict_test,CN_pDC_test$CN)

gc()
rm(xgb_pDC_CN)
















###__________________________________________________________________________________________###
###____________________________________________Mono_classical________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_classical <- read.delim("~/Analysis/data/protein/Mono_classical_CN_log10.csv",sep=";",dec = ",")
CN_Mono_classical <- subset(CN_Mono_classical,CN_Mono_classical$CN>0)
dim(CN_Mono_classical)


TPM_Mono_classical <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_classical_TPM_log10.csv",sep=";",dec=",")
TPM_Mono_classical$TPM <- 10^TPM_Mono_classical$TPM
TPM_Mono_classical <- merge(TPM_Mono_classical,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Mono_classical <- merge(TPM_Mono_classical,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Mono_classical <- ddply(TPM_Mono_classical,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Mono_classical$TPM <- log10(TPM_Mono_classical$TPM)
#dim(TPM_Mono_classical[duplicated(TPM_Mono_classical$uniprotswissprot),])



CN_Mono_classical <- merge(CN_Mono_classical,TPM_Mono_classical,by.x="ID",by.y="uniprotswissprot")
CN_Mono_classical$ensembl_gene_id <- NULL
CN_Mono_classical <- merge(CN_Mono_classical,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_classical[is.na(CN_Mono_classical)]= 0
print(dim(CN_Mono_classical)) # 7984 7126

CN_Mono_classical_dups <- CN_Mono_classical[duplicated(CN_Mono_classical$ID),]

rownames(CN_Mono_classical) <- CN_Mono_classical$ID
CN_Mono_classical$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Mono_classical[is.na(CN_Mono_classical)]=0
CN_Mono_classical <- CN_Mono_classical[, -nearZeroVar(CN_Mono_classical, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Mono_classical)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Mono_classical), 0.8 * nrow(CN_Mono_classical))
CN_Mono_classical_train <- CN_Mono_classical[train_row,]
CN_Mono_classical_test <- CN_Mono_classical[-train_row,]

dim(CN_Mono_classical_train)
dim(CN_Mono_classical_test)

write.table(CN_Mono_classical_train,"~/Analysis/split_dataset/CN_Mono_classical_with_RNA_train",sep=";")
write.table(CN_Mono_classical_test,"~/Analysis/split_dataset/CN_Mono_classical_with_RNA_test",sep=";")

# CN_Mono_classical_train <- read.delim("~/Analysis/split_dataset/CN_Mono_classical_with_RNA_train",sep=";")
# CN_Mono_classical_test <- read.delim("~/Analysis/split_dataset/CN_Mono_classical_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Mono_classical_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_Mono_classical_CN <- train(CN~.,
                         data=CN_Mono_classical_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_Mono_classical_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_classical_CN, "~/Analysis/models/protein/xgb_CN_Mono_classical_with_RNA_02_06_2022.RDS")
# xgb_Mono_classical_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_classical_with_RNA_02_06_2022.RDS")
#
# print(xgb_Mono_classical_CN)

xgb_Mono_classical_CN_predict_test <- predict(xgb_Mono_classical_CN,CN_Mono_classical_test)
xgb_Mono_classical_CN_predict_test <- data.frame(xgb_Mono_classical_CN_predict_test)
# print(lm(xgb_Mono_classical_CN_predict_test$xgb_Mono_classical_CN_predict_test~CN_Mono_classical_test$CN))
print(cor(xgb_Mono_classical_CN_predict_test$xgb_Mono_classical_CN_predict_test,CN_Mono_classical_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Mono_classical_with_RNA_train$TPM,CN_Mono_classical_with_RNA_train$CN)
print(cor(CN_Mono_classical_train$TPM,CN_Mono_classical_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Mono_classical_CN_predict_test$xgb_Mono_classical_CN_predict_test,CN_Mono_classical_test$CN)

gc()
rm(xgb_Mono_classical_CN)
















###__________________________________________________________________________________________###
###____________________________________________Mono_non_classical________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_non_classical <- read.delim("~/Analysis/data/protein/Mono_non_classical_CN_log10.csv",sep=";",dec = ",")
CN_Mono_non_classical <- subset(CN_Mono_non_classical,CN_Mono_non_classical$CN>0)
dim(CN_Mono_non_classical)


TPM_Mono_non_classical <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_non_classical_TPM_log10.csv",sep=";",dec=",")
TPM_Mono_non_classical$TPM <- 10^TPM_Mono_non_classical$TPM
TPM_Mono_non_classical <- merge(TPM_Mono_non_classical,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Mono_non_classical <- merge(TPM_Mono_non_classical,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Mono_non_classical <- ddply(TPM_Mono_non_classical,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Mono_non_classical$TPM <- log10(TPM_Mono_non_classical$TPM)
#dim(TPM_Mono_non_classical[duplicated(TPM_Mono_non_classical$uniprotswissprot),])



CN_Mono_non_classical <- merge(CN_Mono_non_classical,TPM_Mono_non_classical,by.x="ID",by.y="uniprotswissprot")
CN_Mono_non_classical$ensembl_gene_id <- NULL
CN_Mono_non_classical <- merge(CN_Mono_non_classical,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_non_classical[is.na(CN_Mono_non_classical)]= 0
print(dim(CN_Mono_non_classical)) # 7984 7126

CN_Mono_non_classical_dups <- CN_Mono_non_classical[duplicated(CN_Mono_non_classical$ID),]

rownames(CN_Mono_non_classical) <- CN_Mono_non_classical$ID
CN_Mono_non_classical$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Mono_non_classical[is.na(CN_Mono_non_classical)]=0
CN_Mono_non_classical <- CN_Mono_non_classical[, -nearZeroVar(CN_Mono_non_classical, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Mono_non_classical)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Mono_non_classical), 0.8 * nrow(CN_Mono_non_classical))
CN_Mono_non_classical_train <- CN_Mono_non_classical[train_row,]
CN_Mono_non_classical_test <- CN_Mono_non_classical[-train_row,]

dim(CN_Mono_non_classical_train)
dim(CN_Mono_non_classical_test)

write.table(CN_Mono_non_classical_train,"~/Analysis/split_dataset/CN_Mono_non_classical_with_RNA_train",sep=";")
write.table(CN_Mono_non_classical_test,"~/Analysis/split_dataset/CN_Mono_non_classical_with_RNA_test",sep=";")

# CN_Mono_non_classical_train <- read.delim("~/Analysis/split_dataset/CN_Mono_non_classical_with_RNA_train",sep=";")
# CN_Mono_non_classical_test <- read.delim("~/Analysis/split_dataset/CN_Mono_non_classical_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Mono_non_classical_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_Mono_non_classical_CN <- train(CN~.,
                         data=CN_Mono_non_classical_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_Mono_non_classical_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_non_classical_CN, "~/Analysis/models/protein/xgb_CN_Mono_non_classical_with_RNA_02_06_2022.RDS")
# xgb_Mono_non_classical_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_non_classical_with_RNA_02_06_2022.RDS")
#
# print(xgb_Mono_non_classical_CN)

xgb_Mono_non_classical_CN_predict_test <- predict(xgb_Mono_non_classical_CN,CN_Mono_non_classical_test)
xgb_Mono_non_classical_CN_predict_test <- data.frame(xgb_Mono_non_classical_CN_predict_test)
# print(lm(xgb_Mono_non_classical_CN_predict_test$xgb_Mono_non_classical_CN_predict_test~CN_Mono_non_classical_test$CN))
print(cor(xgb_Mono_non_classical_CN_predict_test$xgb_Mono_non_classical_CN_predict_test,CN_Mono_non_classical_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Mono_non_classical_with_RNA_train$TPM,CN_Mono_non_classical_with_RNA_train$CN)
print(cor(CN_Mono_non_classical_train$TPM,CN_Mono_non_classical_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Mono_non_classical_CN_predict_test$xgb_Mono_non_classical_CN_predict_test,CN_Mono_non_classical_test$CN)

gc()
rm(xgb_Mono_non_classical_CN)















###__________________________________________________________________________________________###
###____________________________________________Mono_inter________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_inter <- read.delim("~/Analysis/data/protein/Mono_intermediate_CN_log10.csv",sep=";",dec = ",")
CN_Mono_inter <- subset(CN_Mono_inter,CN_Mono_inter$CN>0)
dim(CN_Mono_inter)


TPM_Mono_inter <- read.delim("~/Analysis/data/immune_cells_RNA/Mono_inter_TPM_log10.csv",sep=";",dec=",")
TPM_Mono_inter$TPM <- 10^TPM_Mono_inter$TPM
TPM_Mono_inter <- merge(TPM_Mono_inter,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Mono_inter <- merge(TPM_Mono_inter,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Mono_inter <- ddply(TPM_Mono_inter,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Mono_inter$TPM <- log10(TPM_Mono_inter$TPM)
#dim(TPM_Mono_inter[duplicated(TPM_Mono_inter$uniprotswissprot),])



CN_Mono_inter <- merge(CN_Mono_inter,TPM_Mono_inter,by.x="ID",by.y="uniprotswissprot")
CN_Mono_inter$ensembl_gene_id <- NULL
CN_Mono_inter <- merge(CN_Mono_inter,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_inter[is.na(CN_Mono_inter)]= 0
print(dim(CN_Mono_inter)) # 7984 7126

CN_Mono_inter_dups <- CN_Mono_inter[duplicated(CN_Mono_inter$ID),]

rownames(CN_Mono_inter) <- CN_Mono_inter$ID
CN_Mono_inter$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Mono_inter[is.na(CN_Mono_inter)]=0
CN_Mono_inter <- CN_Mono_inter[, -nearZeroVar(CN_Mono_inter, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Mono_inter)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Mono_inter), 0.8 * nrow(CN_Mono_inter))
CN_Mono_inter_train <- CN_Mono_inter[train_row,]
CN_Mono_inter_test <- CN_Mono_inter[-train_row,]

dim(CN_Mono_inter_train)
dim(CN_Mono_inter_test)

write.table(CN_Mono_inter_train,"~/Analysis/split_dataset/CN_Mono_inter_with_RNA_train",sep=";")
write.table(CN_Mono_inter_test,"~/Analysis/split_dataset/CN_Mono_inter_with_RNA_test",sep=";")

CN_Mono_inter_train <- read.delim("~/Analysis/split_dataset/CN_Mono_inter_with_RNA_train",sep=";")
CN_Mono_inter_test <- read.delim("~/Analysis/split_dataset/CN_Mono_inter_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Mono_inter_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                   max_depth = 6,
                                   colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                   ## The values below are default values in the sklearn-api.
                                   eta = 0.05,
                                   gamma=0,
                                   min_child_weight = 0.9,
                                   subsample = 1)

start_time <- Sys.time()
xgb_Mono_inter_CN <- train(CN~.,
                         data=CN_Mono_inter_train,
                         method="xgbTree",
                         trControl=control,
                         #metric="Rsquared",
                         tuneGrid= xgbGrid_Mono_inter_CN,
                         na.action = na.omit,
                         nthread=90,
                         verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_inter_CN, "~/Analysis/models/protein/xgb_CN_Mono_inter_with_RNA_02_06_2022.RDS")
# xgb_Mono_inter_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_inter_with_RNA_02_06_2022.RDS")
#
# print(xgb_Mono_inter_CN)

xgb_Mono_inter_CN_predict_test <- predict(xgb_Mono_inter_CN,CN_Mono_inter_test)
xgb_Mono_inter_CN_predict_test <- data.frame(xgb_Mono_inter_CN_predict_test)
# print(lm(xgb_Mono_inter_CN_predict_test$xgb_Mono_inter_CN_predict_test~CN_Mono_inter_test$CN))
print(cor(xgb_Mono_inter_CN_predict_test$xgb_Mono_inter_CN_predict_test,CN_Mono_inter_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Mono_inter_with_RNA_train$TPM,CN_Mono_inter_with_RNA_train$CN)
print(cor(CN_Mono_inter_train$TPM,CN_Mono_inter_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Mono_inter_CN_predict_test$xgb_Mono_inter_CN_predict_test,CN_Mono_inter_test$CN)

gc()
rm(xgb_Mono_inter_CN)
















###__________________________________________________________________________________________###
###____________________________________________Tregs________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_nTregs <- read.delim("~/Analysis/data/protein/nTregs_CN_log10.csv",sep=";",dec = ",")
CN_mTregs <- read.delim("~/Analysis/data/protein/mTregs_CN_log10.csv",sep=";",dec = ",")
CN_Tregs <- merge(CN_nTregs,CN_nTregs,by="ID", all=T)

CN_Tregs$CN <- (10^CN_Tregs$CN.x + 10^CN_Tregs$CN.y)/2
CN_Tregs$CN <- log10(CN_Tregs$CN)
CN_Tregs$CN.x <- NULL
CN_Tregs$CN.y <- NULL

CN_Tregs <- subset(CN_Tregs,CN_Tregs$CN>0)
dim(CN_Tregs)


TPM_Tregs <- read.delim("~/Analysis/data/immune_cells_RNA/Tregs_TPM_log10.csv",sep=";",dec=",")


TPM_Tregs <- merge(TPM_Tregs,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_Tregs <- merge(TPM_Tregs,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_Tregs <- ddply(TPM_Tregs,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_Tregs$TPM <- log10(TPM_Tregs$TPM)
#dim(TPM_Tregs[duplicated(TPM_Tregs$uniprotswissprot),])



CN_Tregs <- merge(CN_Tregs,TPM_Tregs,by.x="ID",by.y="uniprotswissprot")
CN_Tregs$ensembl_gene_id <- NULL
CN_Tregs <- merge(CN_Tregs,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Tregs[is.na(CN_Tregs)]= 0
print(dim(CN_Tregs)) # 7984 7126

CN_Tregs_dups <- CN_Tregs[duplicated(CN_Tregs$ID),]

rownames(CN_Tregs) <- CN_Tregs$ID
CN_Tregs$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Tregs[is.na(CN_Tregs)]=0
CN_Tregs <- CN_Tregs[, -nearZeroVar(CN_Tregs, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Tregs)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Tregs), 0.8 * nrow(CN_Tregs))
CN_Tregs_train <- CN_Tregs[train_row,]
CN_Tregs_test <- CN_Tregs[-train_row,]

dim(CN_Tregs_train)
dim(CN_Tregs_test)

write.table(CN_Tregs_train,"~/Analysis/split_dataset/CN_Tregs_with_RNA_train",sep=";")
write.table(CN_Tregs_test,"~/Analysis/split_dataset/CN_Tregs_with_RNA_test",sep=";")

# CN_Tregs_train <- read.delim("~/Analysis/split_dataset/CN_Tregs_with_RNA_train",sep=";")
# CN_Tregs_test <- read.delim("~/Analysis/split_dataset/CN_Tregs_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Tregs_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                max_depth = 6,
                                colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                ## The values below are default values in the sklearn-api.
                                eta = 0.05,
                                gamma=0,
                                min_child_weight = 0.9,
                                subsample = 1)

start_time <- Sys.time()
xgb_Tregs_CN <- train(CN~.,
                      data=CN_Tregs_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid_Tregs_CN,
                      na.action = na.omit,
                      nthread=90,
                      verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Tregs_CN, "~/Analysis/models/protein/xgb_CN_Tregs_with_RNA_02_06_2022.RDS")
# xgb_Tregs_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Tregs_with_RNA_02_06_2022.RDS")
#
# print(xgb_Tregs_CN)

xgb_Tregs_CN_predict_test <- predict(xgb_Tregs_CN,CN_Tregs_test)
xgb_Tregs_CN_predict_test <- data.frame(xgb_Tregs_CN_predict_test)
# print(lm(xgb_Tregs_CN_predict_test$xgb_Tregs_CN_predict_test~CN_Tregs_test$CN))
print(cor(xgb_Tregs_CN_predict_test$xgb_Tregs_CN_predict_test,CN_Tregs_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_Tregs_with_RNA_train$TPM,CN_Tregs_with_RNA_train$CN)
print(cor(CN_Tregs_train$TPM,CN_Tregs_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_Tregs_CN_predict_test$xgb_Tregs_CN_predict_test,CN_Tregs_test$CN)

gc()
rm(xgb_Tregs_CN)

















###__________________________________________________________________________________________###
###____________________________________________NK________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_NK_bright <- read.delim("~/Analysis/data/protein/NK_bright_CN_log10.csv",sep=";",dec = ",")
CN_NK_dim <- read.delim("~/Analysis/data/protein/NK_dim_CN_log10.csv",sep=";",dec = ",")
CN_NK <- merge(CN_NK_bright,CN_NK_dim,by="ID", all=T)

CN_NK$CN <- (10^CN_NK$CN.x + 10^CN_NK$CN.y)/2
CN_NK$CN <- log10(CN_NK$CN)
CN_NK$CN.x <- NULL
CN_NK$CN.y <- NULL

CN_NK <- subset(CN_NK,CN_NK$CN>0)
dim(CN_NK)


TPM_NK <- read.delim("~/Analysis/data/immune_cells_RNA/NK_TPM_log10.csv",sep=";",dec=",")


TPM_NK <- merge(TPM_NK,tx2gene, by.x="ID", by.y="ensembl_transcript_id")


registerDoMC(10)

TPM_NK <- merge(TPM_NK,gene2PtID,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
TPM_NK <- ddply(TPM_NK,"uniprotswissprot", numcolwise(sum), .parallel = T, .progress = T)
TPM_NK$TPM <- log10(TPM_NK$TPM)
#dim(TPM_NK[duplicated(TPM_NK$uniprotswissprot),])



CN_NK <- merge(CN_NK,TPM_NK,by.x="ID",by.y="uniprotswissprot")
CN_NK$ensembl_gene_id <- NULL
CN_NK <- merge(CN_NK,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_NK[is.na(CN_NK)]= 0
print(dim(CN_NK)) # 7984 7126

CN_NK_dups <- CN_NK[duplicated(CN_NK$ID),]

rownames(CN_NK) <- CN_NK$ID
CN_NK$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_NK[is.na(CN_NK)]=0
CN_NK <- CN_NK[, -nearZeroVar(CN_NK, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_NK)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_NK), 0.8 * nrow(CN_NK))
CN_NK_train <- CN_NK[train_row,]
CN_NK_test <- CN_NK[-train_row,]

dim(CN_NK_train)
dim(CN_NK_test)

write.table(CN_NK_train,"~/Analysis/split_dataset/CN_NK_with_RNA_train",sep=";")
write.table(CN_NK_test,"~/Analysis/split_dataset/CN_NK_with_RNA_test",sep=";")

# CN_NK_train <- read.delim("~/Analysis/split_dataset/CN_NK_with_RNA_train",sep=";")
# CN_NK_test <- read.delim("~/Analysis/split_dataset/CN_NK_with_RNA_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling CD8 TN CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_NK_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                max_depth = 6,
                                colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                ## The values below are default values in the sklearn-api.
                                eta = 0.05,
                                gamma=0,
                                min_child_weight = 0.9,
                                subsample = 1)

start_time <- Sys.time()
xgb_NK_CN <- train(CN~.,
                      data=CN_NK_train,
                      method="xgbTree",
                      trControl=control,
                      #metric="Rsquared",
                      tuneGrid= xgbGrid_NK_CN,
                      na.action = na.omit,
                      nthread=90,
                      verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_NK_CN, "~/Analysis/models/protein/xgb_CN_NK_with_RNA_02_06_2022.RDS")
# xgb_NK_CN <- readRDS("~/Analysis/models/protein/xgb_CN_NK_with_RNA_02_06_2022.RDS")
#
# print(xgb_NK_CN)

xgb_NK_CN_predict_test <- predict(xgb_NK_CN,CN_NK_test)
xgb_NK_CN_predict_test <- data.frame(xgb_NK_CN_predict_test)
# print(lm(xgb_NK_CN_predict_test$xgb_NK_CN_predict_test~CN_NK_test$CN))
print(cor(xgb_NK_CN_predict_test$xgb_NK_CN_predict_test,CN_NK_test$CN, method = "pearson", use = "complete.obs")^2)

# plot(CN_NK_with_RNA_train$TPM,CN_NK_with_RNA_train$CN)
print(cor(CN_NK_train$TPM,CN_NK_train$CN, method = "pearson", use = "complete.obs")^2)

# plot(xgb_NK_CN_predict_test$xgb_NK_CN_predict_test,CN_NK_test$CN)

gc()
rm(xgb_NK_CN)








