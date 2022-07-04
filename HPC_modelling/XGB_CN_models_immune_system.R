library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)

setwd("~/Analysis/")


###__________________________________________________________________________________________###
###____________________________________________ex vivo________________________________________###
###__________________________________________________________________________________________###


## parameters ##
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")



###__________________________________________________________________________________________###
###____________________________________________B naive________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_B_naive <- read.delim("~/Analysis/data/protein/B_naive_CN_log10.csv",sep=";",dec = ",")
CN_B_naive <- subset(CN_B_naive,CN_B_naive$CN>0)
dim(CN_B_naive)

CN_B_naive <- merge(CN_B_naive,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_B_naive[is.na(CN_B_naive)]= 0
print(dim(CN_B_naive)) # 7984 7126


rownames(CN_B_naive) <- CN_B_naive$ID
CN_B_naive$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_B_naive[is.na(CN_B_naive)]=0
CN_B_naive <- CN_B_naive[, -nearZeroVar(CN_B_naive, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_B_naive)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_B_naive), 0.8 * nrow(CN_B_naive))
CN_B_naive_train <- CN_B_naive[train_row,]
CN_B_naive_test <- CN_B_naive[-train_row,]

dim(CN_B_naive_train)
# dim(CN_B_naive_test)

# write.table(CN_B_naive_train,"~/Analysis/split_dataset/CN_B_naive_train",sep=";")
# write.table(CN_B_naive_test,"~/Analysis/split_dataset/CN_B_naive_test",sep=";")

CN_B_naive_train <- read.delim("~/Analysis/split_dataset/CN_B_naive_train",sep=";")
CN_B_naive_test <- read.delim("~/Analysis/split_dataset/CN_B_naive_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling B_naive CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_B_naive_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                              max_depth = 6,
                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                              ## The values below are default values in the sklearn-api.
                              eta = 0.05,
                              gamma=0,
                              min_child_weight = 0.9,
                              subsample = 1)

start_time <- Sys.time()
xgb_B_naive_CN <- train(CN~.,
                        data=CN_B_naive_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_B_naive_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_B_naive_CN, "~/Analysis/models/protein/xgb_CN_B_naive_02_05_2022.RDS")
# xgb_B_naive_CN <- readRDS("~/Analysis/models/protein/xgb_CN_B_naive_02_05_2022.RDS")
#
# print(xgb_B_naive_CN)

xgb_B_naive_CN_predict_test <- predict(xgb_B_naive_CN,CN_B_naive_test)
xgb_B_naive_CN_predict_test <- data.frame(xgb_B_naive_CN_predict_test)
print(lm(xgb_B_naive_CN_predict_test$xgb_B_naive_CN_predict_test~CN_B_naive_test$CN))
print(cor(xgb_B_naive_CN_predict_test$xgb_B_naive_CN_predict_test,CN_B_naive_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_B_naive_CN)
















###__________________________________________________________________________________________###
###____________________________________________B mem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_B_mem <- read.delim("~/Analysis/data/protein/B_mem_CN_log10.csv",sep=";",dec = ",")
CN_B_mem <- subset(CN_B_mem,CN_B_mem$CN>0)
dim(CN_B_mem)

CN_B_mem <- merge(CN_B_mem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_B_mem[is.na(CN_B_mem)]= 0
print(dim(CN_B_mem)) # 7984 7126


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
# dim(CN_B_mem_test)

write.table(CN_B_mem_train,"~/Analysis/split_dataset/CN_B_mem_train",sep=";")
write.table(CN_B_mem_test,"~/Analysis/split_dataset/CN_B_mem_test",sep=";")

CN_B_mem_train <- read.delim("~/Analysis/split_dataset/CN_B_mem_train",sep=";")
CN_B_mem_test <- read.delim("~/Analysis/split_dataset/CN_B_mem_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling B_mem CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_B_mem_CN, "~/Analysis/models/protein/xgb_CN_B_mem_02_05_2022.RDS")
# xgb_B_mem_CN <- readRDS("~/Analysis/models/protein/xgb_CN_B_mem_02_05_2022.RDS")
#
# print(xgb_B_mem_CN)

xgb_B_mem_CN_predict_test <- predict(xgb_B_mem_CN,CN_B_mem_test)
xgb_B_mem_CN_predict_test <- data.frame(xgb_B_mem_CN_predict_test)
print(lm(xgb_B_mem_CN_predict_test$xgb_B_mem_CN_predict_test~CN_B_mem_test$CN))
print(cor(xgb_B_mem_CN_predict_test$xgb_B_mem_CN_predict_test,CN_B_mem_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_B_mem_CN)


















###__________________________________________________________________________________________###
###____________________________________________B plasma________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_B_plasma <- read.delim("~/Analysis/data/protein/B_plasma_CN_log10.csv",sep=";",dec = ",")
CN_B_plasma <- subset(CN_B_plasma,CN_B_plasma$CN>0)
dim(CN_B_plasma)

CN_B_plasma <- merge(CN_B_plasma,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_B_plasma[is.na(CN_B_plasma)]= 0
print(dim(CN_B_plasma)) # 7984 7126


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
# dim(CN_B_plasma_test)

write.table(CN_B_plasma_train,"~/Analysis/split_dataset/CN_B_plasma_train",sep=";")
write.table(CN_B_plasma_test,"~/Analysis/split_dataset/CN_B_plasma_test",sep=";")

CN_B_plasma_train <- read.delim("~/Analysis/split_dataset/CN_B_plasma_train",sep=";")
CN_B_plasma_test <- read.delim("~/Analysis/split_dataset/CN_B_plasma_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling B_plasma CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_B_plasma_CN, "~/Analysis/models/protein/xgb_CN_B_plasma_02_05_2022.RDS")
# xgb_B_plasma_CN <- readRDS("~/Analysis/models/protein/xgb_CN_B_plasma_02_05_2022.RDS")
#
# print(xgb_B_plasma_CN)

xgb_B_plasma_CN_predict_test <- predict(xgb_B_plasma_CN,CN_B_plasma_test)
xgb_B_plasma_CN_predict_test <- data.frame(xgb_B_plasma_CN_predict_test)
print(lm(xgb_B_plasma_CN_predict_test$xgb_B_plasma_CN_predict_test~CN_B_plasma_test$CN))
print(cor(xgb_B_plasma_CN_predict_test$xgb_B_plasma_CN_predict_test,CN_B_plasma_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_B_plasma_CN)












###__________________________________________________________________________________________###
###____________________________________________Basophils________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Basophils <- read.delim("~/Analysis/data/protein/Basophils_CN_log10.csv",sep=";",dec = ",")
CN_Basophils <- subset(CN_Basophils,CN_Basophils$CN>0)
dim(CN_Basophils)

CN_Basophils <- merge(CN_Basophils,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Basophils[is.na(CN_Basophils)]= 0
print(dim(CN_Basophils)) # 7984 7126


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
# dim(CN_Basophils_test)

write.table(CN_Basophils_train,"~/Analysis/split_dataset/CN_Basophils_train",sep=";")
write.table(CN_Basophils_test,"~/Analysis/split_dataset/CN_Basophils_test",sep=";")

CN_Basophils_train <- read.delim("~/Analysis/split_dataset/CN_Basophils_train",sep=";")
CN_Basophils_test <- read.delim("~/Analysis/split_dataset/CN_Basophils_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Basophils CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Basophils_CN, "~/Analysis/models/protein/xgb_CN_Basophils_02_05_2022.RDS")
# xgb_Basophils_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Basophils_02_05_2022.RDS")
#
# print(xgb_Basophils_CN)

xgb_Basophils_CN_predict_test <- predict(xgb_Basophils_CN,CN_Basophils_test)
xgb_Basophils_CN_predict_test <- data.frame(xgb_Basophils_CN_predict_test)
print(lm(xgb_Basophils_CN_predict_test$xgb_Basophils_CN_predict_test~CN_Basophils_test$CN))
print(cor(xgb_Basophils_CN_predict_test$xgb_Basophils_CN_predict_test,CN_Basophils_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Basophils_CN)












###__________________________________________________________________________________________###
###________________________________________ Eosinophils _____________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Eosinophils <- read.delim("~/Analysis/data/protein/Eosinophils_CN_log10.csv",sep=";",dec = ",")
CN_Eosinophils <- subset(CN_Eosinophils,CN_Eosinophils$CN>0)
dim(CN_Eosinophils)

CN_Eosinophils <- merge(CN_Eosinophils,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Eosinophils[is.na(CN_Eosinophils)]= 0
print(dim(CN_Eosinophils)) # 7984 7126


rownames(CN_Eosinophils) <- CN_Eosinophils$ID
CN_Eosinophils$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Eosinophils[is.na(CN_Eosinophils)]=0
CN_Eosinophils <- CN_Eosinophils[, -nearZeroVar(CN_Eosinophils, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Eosinophils)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Eosinophils), 0.8 * nrow(CN_Eosinophils))
CN_Eosinophils_train <- CN_Eosinophils[train_row,]
CN_Eosinophils_test <- CN_Eosinophils[-train_row,]

dim(CN_Eosinophils_train)
# dim(CN_Eosinophils_test)

write.table(CN_Eosinophils_train,"~/Analysis/split_dataset/CN_Eosinophils_train",sep=";")
write.table(CN_Eosinophils_test,"~/Analysis/split_dataset/CN_Eosinophils_test",sep=";")

CN_Eosinophils_train <- read.delim("~/Analysis/split_dataset/CN_Eosinophils_train",sep=";")
CN_Eosinophils_test <- read.delim("~/Analysis/split_dataset/CN_Eosinophils_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Eosinophils CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Eosinophils_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_Eosinophils_CN <- train(CN~.,
                        data=CN_Eosinophils_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_Eosinophils_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Eosinophils_CN, "~/Analysis/models/protein/xgb_CN_Eosinophils_02_05_2022.RDS")
# xgb_Eosinophils_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Eosinophils_02_05_2022.RDS")
#
# print(xgb_Eosinophils_CN)

xgb_Eosinophils_CN_predict_test <- predict(xgb_Eosinophils_CN,CN_Eosinophils_test)
xgb_Eosinophils_CN_predict_test <- data.frame(xgb_Eosinophils_CN_predict_test)
print(lm(xgb_Eosinophils_CN_predict_test$xgb_Eosinophils_CN_predict_test~CN_Eosinophils_test$CN))
print(cor(xgb_Eosinophils_CN_predict_test$xgb_Eosinophils_CN_predict_test,CN_Eosinophils_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Eosinophils_CN)











###__________________________________________________________________________________________###
###____________________________________________ mDC  ________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_mDC <- read.delim("~/Analysis/data/protein/mDC_CN_log10.csv",sep=";",dec = ",")
CN_mDC <- subset(CN_mDC,CN_mDC$CN>0)
dim(CN_mDC)

CN_mDC <- merge(CN_mDC,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_mDC[is.na(CN_mDC)]= 0
print(dim(CN_mDC)) # 7984 7126


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
# dim(CN_mDC_test)

write.table(CN_mDC_train,"~/Analysis/split_dataset/CN_mDC_train",sep=";")
write.table(CN_mDC_test,"~/Analysis/split_dataset/CN_mDC_test",sep=";")

CN_mDC_train <- read.delim("~/Analysis/split_dataset/CN_mDC_train",sep=";")
CN_mDC_test <- read.delim("~/Analysis/split_dataset/CN_mDC_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling mDC CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_mDC_CN, "~/Analysis/models/protein/xgb_CN_mDC_02_05_2022.RDS")
# xgb_mDC_CN <- readRDS("~/Analysis/models/protein/xgb_CN_mDC_02_05_2022.RDS")
#
# print(xgb_mDC_CN)

xgb_mDC_CN_predict_test <- predict(xgb_mDC_CN,CN_mDC_test)
xgb_mDC_CN_predict_test <- data.frame(xgb_mDC_CN_predict_test)
print(lm(xgb_mDC_CN_predict_test$xgb_mDC_CN_predict_test~CN_mDC_test$CN))
print(cor(xgb_mDC_CN_predict_test$xgb_mDC_CN_predict_test,CN_mDC_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_mDC_CN)













###__________________________________________________________________________________________###
###________________________________________ Mono Classical __________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_classical <- read.delim("~/Analysis/data/protein/Mono_classical_CN_log10.csv",sep=";",dec = ",")
CN_Mono_classical <- subset(CN_Mono_classical,CN_Mono_classical$CN>0)
dim(CN_Mono_classical)

CN_Mono_classical <- merge(CN_Mono_classical,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_classical[is.na(CN_Mono_classical)]= 0
print(dim(CN_Mono_classical)) # 7984 7126


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
# dim(CN_Mono_classical_test)

write.table(CN_Mono_classical_train,"~/Analysis/split_dataset/CN_Mono_classical_train",sep=";")
write.table(CN_Mono_classical_test,"~/Analysis/split_dataset/CN_Mono_classical_test",sep=";")

CN_Mono_classical_train <- read.delim("~/Analysis/split_dataset/CN_Mono_classical_train",sep=";")
CN_Mono_classical_test <- read.delim("~/Analysis/split_dataset/CN_Mono_classical_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Mono_classical CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_classical_CN, "~/Analysis/models/protein/xgb_CN_Mono_classical_02_05_2022.RDS")
# xgb_Mono_classical_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_classical_02_05_2022.RDS")
#
# print(xgb_Mono_classical_CN)

xgb_Mono_classical_CN_predict_test <- predict(xgb_Mono_classical_CN,CN_Mono_classical_test)
xgb_Mono_classical_CN_predict_test <- data.frame(xgb_Mono_classical_CN_predict_test)
print(lm(xgb_Mono_classical_CN_predict_test$xgb_Mono_classical_CN_predict_test~CN_Mono_classical_test$CN))
print(cor(xgb_Mono_classical_CN_predict_test$xgb_Mono_classical_CN_predict_test,CN_Mono_classical_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Mono_classical_CN)













###__________________________________________________________________________________________###
###____________________________________ Mono intermediate ___________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_intermediate <- read.delim("~/Analysis/data/protein/Mono_intermediate_CN_log10.csv",sep=";",dec = ",")
CN_Mono_intermediate <- subset(CN_Mono_intermediate,CN_Mono_intermediate$CN>0)
dim(CN_Mono_intermediate)

CN_Mono_intermediate <- merge(CN_Mono_intermediate,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_intermediate[is.na(CN_Mono_intermediate)]= 0
print(dim(CN_Mono_intermediate)) # 7984 7126


rownames(CN_Mono_intermediate) <- CN_Mono_intermediate$ID
CN_Mono_intermediate$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_Mono_intermediate[is.na(CN_Mono_intermediate)]=0
CN_Mono_intermediate <- CN_Mono_intermediate[, -nearZeroVar(CN_Mono_intermediate, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_Mono_intermediate)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_Mono_intermediate), 0.8 * nrow(CN_Mono_intermediate))
CN_Mono_intermediate_train <- CN_Mono_intermediate[train_row,]
CN_Mono_intermediate_test <- CN_Mono_intermediate[-train_row,]

dim(CN_Mono_intermediate_train)
# dim(CN_Mono_intermediate_test)

write.table(CN_Mono_intermediate_train,"~/Analysis/split_dataset/CN_Mono_intermediate_train",sep=";")
write.table(CN_Mono_intermediate_test,"~/Analysis/split_dataset/CN_Mono_intermediate_test",sep=";")

CN_Mono_intermediate_train <- read.delim("~/Analysis/split_dataset/CN_Mono_intermediate_train",sep=";")
CN_Mono_intermediate_test <- read.delim("~/Analysis/split_dataset/CN_Mono_intermediate_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Mono_intermediate CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_Mono_intermediate_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_Mono_intermediate_CN <- train(CN~.,
                        data=CN_Mono_intermediate_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_Mono_intermediate_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_intermediate_CN, "~/Analysis/models/protein/xgb_CN_Mono_intermediate_02_05_2022.RDS")
# xgb_Mono_intermediate_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_intermediate_02_05_2022.RDS")
#
# print(xgb_Mono_intermediate_CN)

xgb_Mono_intermediate_CN_predict_test <- predict(xgb_Mono_intermediate_CN,CN_Mono_intermediate_test)
xgb_Mono_intermediate_CN_predict_test <- data.frame(xgb_Mono_intermediate_CN_predict_test)
print(lm(xgb_Mono_intermediate_CN_predict_test$xgb_Mono_intermediate_CN_predict_test~CN_Mono_intermediate_test$CN))
print(cor(xgb_Mono_intermediate_CN_predict_test$xgb_Mono_intermediate_CN_predict_test,CN_Mono_intermediate_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Mono_intermediate_CN)










###__________________________________________________________________________________________###
###_____________________________________ Mono non-classical _________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Mono_non_classical <- read.delim("~/Analysis/data/protein/Mono_non_classical_CN_log10.csv",sep=";",dec = ",")
CN_Mono_non_classical <- subset(CN_Mono_non_classical,CN_Mono_non_classical$CN>0)
dim(CN_Mono_non_classical)

CN_Mono_non_classical <- merge(CN_Mono_non_classical,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Mono_non_classical[is.na(CN_Mono_non_classical)]= 0
print(dim(CN_Mono_non_classical)) # 7984 7126


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
# dim(CN_Mono_non_classical_test)

write.table(CN_Mono_non_classical_train,"~/Analysis/split_dataset/CN_Mono_non_classical_train",sep=";")
write.table(CN_Mono_non_classical_test,"~/Analysis/split_dataset/CN_Mono_non_classical_test",sep=";")

CN_Mono_non_classical_train <- read.delim("~/Analysis/split_dataset/CN_Mono_non_classical_train",sep=";")
CN_Mono_non_classical_test <- read.delim("~/Analysis/split_dataset/CN_Mono_non_classical_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Mono_non_classical CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Mono_non_classical_CN, "~/Analysis/models/protein/xgb_CN_Mono_non_classical_02_05_2022.RDS")
# xgb_Mono_non_classical_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Mono_non_classical_02_05_2022.RDS")
#
# print(xgb_Mono_non_classical_CN)

xgb_Mono_non_classical_CN_predict_test <- predict(xgb_Mono_non_classical_CN,CN_Mono_non_classical_test)
xgb_Mono_non_classical_CN_predict_test <- data.frame(xgb_Mono_non_classical_CN_predict_test)
print(lm(xgb_Mono_non_classical_CN_predict_test$xgb_Mono_non_classical_CN_predict_test~CN_Mono_non_classical_test$CN))
print(cor(xgb_Mono_non_classical_CN_predict_test$xgb_Mono_non_classical_CN_predict_test,CN_Mono_non_classical_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Mono_non_classical_CN)












###__________________________________________________________________________________________###
###___________________________________________ nTregs _______________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_mTregs <- read.delim("~/Analysis/data/protein/mTregs_CN_log10.csv",sep=";",dec = ",")
CN_mTregs <- subset(CN_mTregs,CN_mTregs$CN>0)
dim(CN_mTregs)

CN_mTregs <- merge(CN_mTregs,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_mTregs[is.na(CN_mTregs)]= 0
print(dim(CN_mTregs)) # 7984 7126


rownames(CN_mTregs) <- CN_mTregs$ID
CN_mTregs$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_mTregs[is.na(CN_mTregs)]=0
CN_mTregs <- CN_mTregs[, -nearZeroVar(CN_mTregs, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_mTregs)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_mTregs), 0.8 * nrow(CN_mTregs))
CN_mTregs_train <- CN_mTregs[train_row,]
CN_mTregs_test <- CN_mTregs[-train_row,]

dim(CN_mTregs_train)
# dim(CN_mTregs_test)

write.table(CN_mTregs_train,"~/Analysis/split_dataset/CN_mTregs_train",sep=";")
write.table(CN_mTregs_test,"~/Analysis/split_dataset/CN_mTregs_test",sep=";")

CN_mTregs_train <- read.delim("~/Analysis/split_dataset/CN_mTregs_train",sep=";")
CN_mTregs_test <- read.delim("~/Analysis/split_dataset/CN_mTregs_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling mTregs CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_mTregs_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_mTregs_CN <- train(CN~.,
                        data=CN_mTregs_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_mTregs_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_mTregs_CN, "~/Analysis/models/protein/xgb_CN_mTregs_02_05_2022.RDS")
# xgb_mTregs_CN <- readRDS("~/Analysis/models/protein/xgb_CN_mTregs_02_05_2022.RDS")
#
# print(xgb_mTregs_CN)

xgb_mTregs_CN_predict_test <- predict(xgb_mTregs_CN,CN_mTregs_test)
xgb_mTregs_CN_predict_test <- data.frame(xgb_mTregs_CN_predict_test)
print(lm(xgb_mTregs_CN_predict_test$xgb_mTregs_CN_predict_test~CN_mTregs_test$CN))
print(cor(xgb_mTregs_CN_predict_test$xgb_mTregs_CN_predict_test,CN_mTregs_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_mTregs_CN)











###__________________________________________________________________________________________###
###___________________________________________ nTregs _______________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_nTregs <- read.delim("~/Analysis/data/protein/nTregs_CN_log10.csv",sep=";",dec = ",")
CN_nTregs <- subset(CN_nTregs,CN_nTregs$CN>0)
dim(CN_nTregs)

CN_nTregs <- merge(CN_nTregs,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_nTregs[is.na(CN_nTregs)]= 0
print(dim(CN_nTregs)) # 7984 7126


rownames(CN_nTregs) <- CN_nTregs$ID
CN_nTregs$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_nTregs[is.na(CN_nTregs)]=0
CN_nTregs <- CN_nTregs[, -nearZeroVar(CN_nTregs, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_nTregs)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_nTregs), 0.8 * nrow(CN_nTregs))
CN_nTregs_train <- CN_nTregs[train_row,]
CN_nTregs_test <- CN_nTregs[-train_row,]

dim(CN_nTregs_train)
# dim(CN_nTregs_test)

write.table(CN_nTregs_train,"~/Analysis/split_dataset/CN_nTregs_train",sep=";")
write.table(CN_nTregs_test,"~/Analysis/split_dataset/CN_nTregs_test",sep=";")

CN_nTregs_train <- read.delim("~/Analysis/split_dataset/CN_nTregs_train",sep=";")
CN_nTregs_test <- read.delim("~/Analysis/split_dataset/CN_nTregs_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling nTregs CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_nTregs_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_nTregs_CN <- train(CN~.,
                        data=CN_nTregs_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_nTregs_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_nTregs_CN, "~/Analysis/models/protein/xgb_CN_nTregs_02_05_2022.RDS")
# xgb_nTregs_CN <- readRDS("~/Analysis/models/protein/xgb_CN_nTregs_02_05_2022.RDS")
#
# print(xgb_nTregs_CN)

xgb_nTregs_CN_predict_test <- predict(xgb_nTregs_CN,CN_nTregs_test)
xgb_nTregs_CN_predict_test <- data.frame(xgb_nTregs_CN_predict_test)
print(lm(xgb_nTregs_CN_predict_test$xgb_nTregs_CN_predict_test~CN_nTregs_test$CN))
print(cor(xgb_nTregs_CN_predict_test$xgb_nTregs_CN_predict_test,CN_nTregs_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_nTregs_CN)










###__________________________________________________________________________________________###
###________________________________________ Neutrophils _____________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_Neutrophils <- read.delim("~/Analysis/data/protein/Neutrophils_CN_log10.csv",sep=";",dec = ",")
CN_Neutrophils <- subset(CN_Neutrophils,CN_Neutrophils$CN>0)
dim(CN_Neutrophils)

CN_Neutrophils <- merge(CN_Neutrophils,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_Neutrophils[is.na(CN_Neutrophils)]= 0
print(dim(CN_Neutrophils)) # 7984 7126


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
# dim(CN_Neutrophils_test)

write.table(CN_Neutrophils_train,"~/Analysis/split_dataset/CN_Neutrophils_train",sep=";")
write.table(CN_Neutrophils_test,"~/Analysis/split_dataset/CN_Neutrophils_test",sep=";")

CN_Neutrophils_train <- read.delim("~/Analysis/split_dataset/CN_Neutrophils_train",sep=";")
CN_Neutrophils_test <- read.delim("~/Analysis/split_dataset/CN_Neutrophils_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling Neutrophils CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_Neutrophils_CN, "~/Analysis/models/protein/xgb_CN_Neutrophils_02_05_2022.RDS")
# xgb_Neutrophils_CN <- readRDS("~/Analysis/models/protein/xgb_CN_Neutrophils_02_05_2022.RDS")
#
# print(xgb_Neutrophils_CN)

xgb_Neutrophils_CN_predict_test <- predict(xgb_Neutrophils_CN,CN_Neutrophils_test)
xgb_Neutrophils_CN_predict_test <- data.frame(xgb_Neutrophils_CN_predict_test)
print(lm(xgb_Neutrophils_CN_predict_test$xgb_Neutrophils_CN_predict_test~CN_Neutrophils_test$CN))
print(cor(xgb_Neutrophils_CN_predict_test$xgb_Neutrophils_CN_predict_test,CN_Neutrophils_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_Neutrophils_CN)









###__________________________________________________________________________________________###
###_________________________________________ NK-bright ______________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_NK_bright <- read.delim("~/Analysis/data/protein/NK_bright_CN_log10.csv",sep=";",dec = ",")
CN_NK_bright <- subset(CN_NK_bright,CN_NK_bright$CN>0)
dim(CN_NK_bright)

CN_NK_bright <- merge(CN_NK_bright,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_NK_bright[is.na(CN_NK_bright)]= 0
print(dim(CN_NK_bright)) # 7984 7126


rownames(CN_NK_bright) <- CN_NK_bright$ID
CN_NK_bright$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_NK_bright[is.na(CN_NK_bright)]=0
CN_NK_bright <- CN_NK_bright[, -nearZeroVar(CN_NK_bright, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_NK_bright)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_NK_bright), 0.8 * nrow(CN_NK_bright))
CN_NK_bright_train <- CN_NK_bright[train_row,]
CN_NK_bright_test <- CN_NK_bright[-train_row,]

dim(CN_NK_bright_train)
# dim(CN_NK_bright_test)

write.table(CN_NK_bright_train,"~/Analysis/split_dataset/CN_NK_bright_train",sep=";")
write.table(CN_NK_bright_test,"~/Analysis/split_dataset/CN_NK_bright_test",sep=";")

CN_NK_bright_train <- read.delim("~/Analysis/split_dataset/CN_NK_bright_train",sep=";")
CN_NK_bright_test <- read.delim("~/Analysis/split_dataset/CN_NK_bright_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling NK_bright CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_NK_bright_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_NK_bright_CN <- train(CN~.,
                        data=CN_NK_bright_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_NK_bright_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_NK_bright_CN, "~/Analysis/models/protein/xgb_CN_NK_bright_02_05_2022.RDS")
# xgb_NK_bright_CN <- readRDS("~/Analysis/models/protein/xgb_CN_NK_bright_02_05_2022.RDS")
#
# print(xgb_NK_bright_CN)

xgb_NK_bright_CN_predict_test <- predict(xgb_NK_bright_CN,CN_NK_bright_test)
xgb_NK_bright_CN_predict_test <- data.frame(xgb_NK_bright_CN_predict_test)
print(lm(xgb_NK_bright_CN_predict_test$xgb_NK_bright_CN_predict_test~CN_NK_bright_test$CN))
print(cor(xgb_NK_bright_CN_predict_test$xgb_NK_bright_CN_predict_test,CN_NK_bright_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_NK_bright_CN)












###__________________________________________________________________________________________###
###___________________________________________ NK-dim _______________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_NK_dim <- read.delim("~/Analysis/data/protein/NK_dim_CN_log10.csv",sep=";",dec = ",")
CN_NK_dim <- subset(CN_NK_dim,CN_NK_dim$CN>0)
dim(CN_NK_dim)

CN_NK_dim <- merge(CN_NK_dim,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_NK_dim[is.na(CN_NK_dim)]= 0
print(dim(CN_NK_dim)) # 7984 7126


rownames(CN_NK_dim) <- CN_NK_dim$ID
CN_NK_dim$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_NK_dim[is.na(CN_NK_dim)]=0
CN_NK_dim <- CN_NK_dim[, -nearZeroVar(CN_NK_dim, allowParallel = T, uniqueCut = 0.1)]
print(dim(CN_NK_dim)) # 7984 2759

gc()






## test / train sets ##
train_row <- sample(1:nrow(CN_NK_dim), 0.8 * nrow(CN_NK_dim))
CN_NK_dim_train <- CN_NK_dim[train_row,]
CN_NK_dim_test <- CN_NK_dim[-train_row,]

dim(CN_NK_dim_train)
# dim(CN_NK_dim_test)

write.table(CN_NK_dim_train,"~/Analysis/split_dataset/CN_NK_dim_train",sep=";")
write.table(CN_NK_dim_test,"~/Analysis/split_dataset/CN_NK_dim_test",sep=";")

CN_NK_dim_train <- read.delim("~/Analysis/split_dataset/CN_NK_dim_train",sep=";")
CN_NK_dim_test <- read.delim("~/Analysis/split_dataset/CN_NK_dim_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling NK_dim CN")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid_NK_dim_CN <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                                  max_depth = 6,
                                  colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                                  ## The values below are default values in the sklearn-api.
                                  eta = 0.05,
                                  gamma=0,
                                  min_child_weight = 0.9,
                                  subsample = 1)

start_time <- Sys.time()
xgb_NK_dim_CN <- train(CN~.,
                        data=CN_NK_dim_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid_NK_dim_CN,
                        na.action = na.omit,
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_NK_dim_CN, "~/Analysis/models/protein/xgb_CN_NK_dim_02_05_2022.RDS")
# xgb_NK_dim_CN <- readRDS("~/Analysis/models/protein/xgb_CN_NK_dim_02_05_2022.RDS")
#
# print(xgb_NK_dim_CN)

xgb_NK_dim_CN_predict_test <- predict(xgb_NK_dim_CN,CN_NK_dim_test)
xgb_NK_dim_CN_predict_test <- data.frame(xgb_NK_dim_CN_predict_test)
print(lm(xgb_NK_dim_CN_predict_test$xgb_NK_dim_CN_predict_test~CN_NK_dim_test$CN))
print(cor(xgb_NK_dim_CN_predict_test$xgb_NK_dim_CN_predict_test,CN_NK_dim_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_NK_dim_CN)











###__________________________________________________________________________________________###
###____________________________________________ pDC _________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_pDC <- read.delim("~/Analysis/data/protein/pDC_CN_log10.csv",sep=";",dec = ",")
CN_pDC <- subset(CN_pDC,CN_pDC$CN>0)
dim(CN_pDC)

CN_pDC <- merge(CN_pDC,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_pDC[is.na(CN_pDC)]= 0
print(dim(CN_pDC)) # 7984 7126


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
# dim(CN_pDC_test)

write.table(CN_pDC_train,"~/Analysis/split_dataset/CN_pDC_train",sep=";")
write.table(CN_pDC_test,"~/Analysis/split_dataset/CN_pDC_test",sep=";")

CN_pDC_train <- read.delim("~/Analysis/split_dataset/CN_pDC_train",sep=";")
CN_pDC_test <- read.delim("~/Analysis/split_dataset/CN_pDC_test",sep=";")
#
#
#
# ## Training the xgb model ##
print("modeling pDC CN")
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
                        nthread=48,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of

saveRDS(xgb_pDC_CN, "~/Analysis/models/protein/xgb_CN_pDC_02_05_2022.RDS")
# xgb_pDC_CN <- readRDS("~/Analysis/models/protein/xgb_CN_pDC_02_05_2022.RDS")
#
# print(xgb_pDC_CN)

xgb_pDC_CN_predict_test <- predict(xgb_pDC_CN,CN_pDC_test)
xgb_pDC_CN_predict_test <- data.frame(xgb_pDC_CN_predict_test)
print(lm(xgb_pDC_CN_predict_test$xgb_pDC_CN_predict_test~CN_pDC_test$CN))
print(cor(xgb_pDC_CN_predict_test$xgb_pDC_CN_predict_test,CN_pDC_test$CN, method = "pearson", use = "complete.obs")^2)
#
#
#
gc()
rm(xgb_pDC_CN)


