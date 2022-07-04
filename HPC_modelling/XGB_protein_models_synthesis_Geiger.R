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
Sequence_parameters_Geiger <- read.delim("~/Analysis/data/synthesis/Sequence_parameters_for_Geiger_data_synthesis_kinectics_14-06-22.csv",sep = ";",dec=",")




# ###__________________________________________________________________________________________###
# ###______________________________________Synthesis_rest______________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# # ## prep data ##
# print("importing and preparing data")
# Synthesis_rest <- read.delim("~/Analysis/data/synthesis/Synthesis_rest_Log2_per_min.csv",sep = ";",dec = ",")
# Synthesis_rest$log2N_per_min <- log10(Synthesis_rest$log2N_per_min)
# 
# Synthesis_rest <- merge(Synthesis_rest,Sequence_parameters_Geiger,by="ID",all.x=F)
# print(dim(Synthesis_rest)) # 6993 6803
# rownames(Synthesis_rest) <- Synthesis_rest$ID
# Synthesis_rest$ID <- NULL
# 
# 
# ## Feature selection ##
# print("Feature selection")
# registerDoMC(4)
# 
# Synthesis_rest[is.na(Synthesis_rest)]=0
# Synthesis_rest <- Synthesis_rest[, -nearZeroVar(Synthesis_rest, allowParallel = T, uniqueCut = 0.1)]
# print(dim(Synthesis_rest)) #
# 
# 
# 
# ## test / train sets ##
# 
# train_row <- sample(1:nrow(Synthesis_rest), 0.8 * nrow(Synthesis_rest))
# Synthesis_rest_train <- Synthesis_rest[train_row,]
# Synthesis_rest_test <- Synthesis_rest[-train_row,]
# 
# print(dim(Synthesis_rest_train))
# print(dim(Synthesis_rest_test))
# 
# write.table(Synthesis_rest_train,"~/Analysis/split_dataset/Synthesis_rest_train.csv",sep=";")
# write.table(Synthesis_rest_test,"~/Analysis/split_dataset/Synthesis_rest_test.csv",sep=";")
# 
# # Synthesis_rest_train <- read.delim("~/Analysis/split_dataset/Synthesis_rest_train.csv",sep=";")
# # Synthesis_rest_test <- read.delim("~/Analysis/split_dataset/Synthesis_rest_test.csv",sep=";")
# 
# 
# 
# 
# gc()
# 
# 
# ## Training the RF model for Tn 0h ##
# print("Modeling synthesis rate rest")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# 
# xgbGrid <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
#                             max_depth = 6,
#                             colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
#                             ## The values below are default values in the sklearn-api.
#                             eta = 0.05,
#                             gamma=0,
#                             min_child_weight = 0.9,
#                             subsample = 1)
# 
# 
# start_time <- Sys.time()
# xgb_Synthesis_rest <- train(log2N_per_min~.,
#                       data=Synthesis_rest_train,
#                       method="xgbTree",
#                       trControl=control,
#                       tuneGrid= xgbGrid,
#                       na.action = na.omit,
#                       nthread=90,
#                       verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time)
# print(xgb_Synthesis_rest)
# 
# 
# saveRDS(xgb_Synthesis_rest, "models/protein/xgb_Synthesis_rest_10k_14_06_2022.RDS")
# # xgb_Synthesis_rest <- readRDS("~/Analysis/models/protein/xgb_Synthesis_rest_10k_14_06_2022.RDS")
# 
# xgb_Synthesis_rest_predict_test <- predict(xgb_Synthesis_rest,Synthesis_rest_test)
# xgb_Synthesis_rest_predict_test <- data.frame(xgb_Synthesis_rest_predict_test)
# 
# plot((xgb_Synthesis_rest_predict_test$xgb_Synthesis_rest_predict_test),(Synthesis_rest_test$log2N_per_min),asp = 1)
# 
# print(cor(log10(xgb_Synthesis_rest_predict_test$xgb_Synthesis_rest_predict_test),log10(Synthesis_rest_test$log2N_per_min),use = "pairwise.complete.obs",method = "pearson"))^2
# 
# gc()
# rm(xgb_Synthesis_rest,Synthesis_rest_train,Synthesis_rest_test)
# 











# ###__________________________________________________________________________________________###
# ###______________________________________Synthesis_6h______________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
print("importing and preparing data")
Synthesis_6h <- read.delim("~/Analysis/data/synthesis/Synthesis_Act_6h_Log2_per_min.csv",sep = ";",dec = ",")
Synthesis_6h$log2N_per_min <- log10(Synthesis_6h$log2N_per_min)

Synthesis_6h <- merge(Synthesis_6h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(Synthesis_6h)) # 6993 6803
rownames(Synthesis_6h) <- Synthesis_6h$ID
Synthesis_6h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

Synthesis_6h[is.na(Synthesis_6h)]=0
Synthesis_6h <- Synthesis_6h[, -nearZeroVar(Synthesis_6h, allowParallel = T, uniqueCut = 0.1)]
print(dim(Synthesis_6h)) #



## test / train sets ##

train_row <- sample(1:nrow(Synthesis_6h), 0.8 * nrow(Synthesis_6h))
Synthesis_6h_train <- Synthesis_6h[train_row,]
Synthesis_6h_test <- Synthesis_6h[-train_row,]

print(dim(Synthesis_6h_train))
print(dim(Synthesis_6h_test))

write.table(Synthesis_6h_train,"~/Analysis/split_dataset/Synthesis_6h_train.csv",sep=";")
write.table(Synthesis_6h_test,"~/Analysis/split_dataset/Synthesis_6h_test.csv",sep=";")

# Synthesis_6h_train <- read.delim("~/Analysis/split_dataset/Synthesis_6h_train.csv",sep=";")
# Synthesis_6h_test <- read.delim("~/Analysis/split_dataset/Synthesis_6h_test.csv",sep=";")




gc()


## Training the RF model for Tn 0h ##
print("Modeling synthesis rate 6h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                       max_depth = 6,
                       colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                       ## The values below are default values in the sklearn-api.
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)


start_time <- Sys.time()
xgb_Synthesis_6h <- train(log2N_per_min~.,
                            data=Synthesis_6h_train,
                            method="xgbTree",
                            trControl=control,
                            tuneGrid= xgbGrid,
                            na.action = na.omit,
                            nthread=90,
                            verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time)
print(xgb_Synthesis_6h)


saveRDS(xgb_Synthesis_6h, "models/protein/xgb_Synthesis_6h_10k_14_06_2022.RDS")
# xgb_Synthesis_6h <- readRDS("~/Analysis/models/protein/xgb_Synthesis_6h_10k_14_06_2022.RDS")

xgb_Synthesis_6h_predict_test <- predict(xgb_Synthesis_6h,Synthesis_6h_test)
xgb_Synthesis_6h_predict_test <- data.frame(xgb_Synthesis_6h_predict_test)

plot((xgb_Synthesis_6h_predict_test$xgb_Synthesis_6h_predict_test),(Synthesis_6h_test$log2N_per_min),asp = 1)

print(cor(log10(xgb_Synthesis_6h_predict_test$xgb_Synthesis_6h_predict_test),log10(Synthesis_6h_test$log2N_per_min),use = "pairwise.complete.obs",method = "pearson"))^2

gc()
rm(xgb_Synthesis_6h,Synthesis_6h_train,Synthesis_6h_test)













# ###__________________________________________________________________________________________###
# ###______________________________________Synthesis_12h______________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
print("importing and preparing data")
Synthesis_12h <- read.delim("~/Analysis/data/synthesis/Synthesis_Act_12h_Log2_per_min.csv",sep = ";",dec = ",")
Synthesis_12h$log2N_per_min <- log10(Synthesis_12h$log2N_per_min)

Synthesis_12h <- merge(Synthesis_12h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(Synthesis_12h)) # 6993 6803
rownames(Synthesis_12h) <- Synthesis_12h$ID
Synthesis_12h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

Synthesis_12h[is.na(Synthesis_12h)]=0
Synthesis_12h <- Synthesis_12h[, -nearZeroVar(Synthesis_12h, allowParallel = T, uniqueCut = 0.1)]
print(dim(Synthesis_12h)) #



## test / train sets ##

train_row <- sample(1:nrow(Synthesis_12h), 0.8 * nrow(Synthesis_12h))
Synthesis_12h_train <- Synthesis_12h[train_row,]
Synthesis_12h_test <- Synthesis_12h[-train_row,]

print(dim(Synthesis_12h_train))
print(dim(Synthesis_12h_test))

write.table(Synthesis_12h_train,"~/Analysis/split_dataset/Synthesis_12h_train.csv",sep=";")
write.table(Synthesis_12h_test,"~/Analysis/split_dataset/Synthesis_12h_test.csv",sep=";")

# Synthesis_12h_train <- read.delim("~/Analysis/split_dataset/Synthesis_12h_train.csv",sep=";")
# Synthesis_12h_test <- read.delim("~/Analysis/split_dataset/Synthesis_12h_test.csv",sep=";")




gc()


## Training the RF model for Tn 0h ##
print("Modeling synthesis rate 12h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                       max_depth = 6,
                       colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                       ## The values below are default values in the sklearn-api.
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)


start_time <- Sys.time()
xgb_Synthesis_12h <- train(log2N_per_min~.,
                            data=Synthesis_12h_train,
                            method="xgbTree",
                            trControl=control,
                            tuneGrid= xgbGrid,
                            na.action = na.omit,
                            nthread=90,
                            verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time)
print(xgb_Synthesis_12h)


saveRDS(xgb_Synthesis_12h, "models/protein/xgb_Synthesis_12h_10k_14_06_2022.RDS")
# xgb_Synthesis_12h <- readRDS("~/Analysis/models/protein/xgb_Synthesis_12h_10k_14_06_2022.RDS")

xgb_Synthesis_12h_predict_test <- predict(xgb_Synthesis_12h,Synthesis_12h_test)
xgb_Synthesis_12h_predict_test <- data.frame(xgb_Synthesis_12h_predict_test)

plot((xgb_Synthesis_12h_predict_test$xgb_Synthesis_12h_predict_test),(Synthesis_12h_test$log2N_per_min),asp = 1)

print(cor(log10(xgb_Synthesis_12h_predict_test$xgb_Synthesis_12h_predict_test),log10(Synthesis_12h_test$log2N_per_min),use = "pairwise.complete.obs",method = "pearson"))^2

gc()
rm(xgb_Synthesis_12h,Synthesis_12h_train,Synthesis_12h_test)















# ###__________________________________________________________________________________________###
# ###______________________________________Synthesis_24h______________________________________###
# ###__________________________________________________________________________________________###
# 
# 
# ## prep data ##
print("importing and preparing data")
Synthesis_24h <- read.delim("~/Analysis/data/synthesis/Synthesis_Act_24h_Log2_per_min.csv",sep = ";",dec = ",")
Synthesis_24h$log2N_per_min <- log10(Synthesis_24h$log2N_per_min)

Synthesis_24h <- merge(Synthesis_24h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(Synthesis_24h)) # 6993 6803
rownames(Synthesis_24h) <- Synthesis_24h$ID
Synthesis_24h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

Synthesis_24h[is.na(Synthesis_24h)]=0
Synthesis_24h <- Synthesis_24h[, -nearZeroVar(Synthesis_24h, allowParallel = T, uniqueCut = 0.1)]
print(dim(Synthesis_24h)) #



## test / train sets ##

train_row <- sample(1:nrow(Synthesis_24h), 0.8 * nrow(Synthesis_24h))
Synthesis_24h_train <- Synthesis_24h[train_row,]
Synthesis_24h_test <- Synthesis_24h[-train_row,]

print(dim(Synthesis_24h_train))
print(dim(Synthesis_24h_test))

write.table(Synthesis_24h_train,"~/Analysis/split_dataset/Synthesis_24h_train.csv",sep=";")
write.table(Synthesis_24h_test,"~/Analysis/split_dataset/Synthesis_24h_test.csv",sep=";")

# Synthesis_24h_train <- read.delim("~/Analysis/split_dataset/Synthesis_24h_train.csv",sep=";")
# Synthesis_24h_test <- read.delim("~/Analysis/split_dataset/Synthesis_24h_test.csv",sep=";")




gc()


## Training the RF model for Tn 0h ##
print("Modeling synthesis rate 24h")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


xgbGrid <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
                       max_depth = 6,
                       colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
                       ## The values below are default values in the sklearn-api.
                       eta = 0.05,
                       gamma=0,
                       min_child_weight = 0.9,
                       subsample = 1)


start_time <- Sys.time()
xgb_Synthesis_24h <- train(log2N_per_min~.,
                            data=Synthesis_24h_train,
                            method="xgbTree",
                            trControl=control,
                            tuneGrid= xgbGrid,
                            na.action = na.omit,
                            nthread=90,
                            verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time)
print(xgb_Synthesis_24h)


saveRDS(xgb_Synthesis_24h, "models/protein/xgb_Synthesis_24h_10k_14_06_2022.RDS")
# xgb_Synthesis_24h <- readRDS("~/Analysis/models/protein/xgb_Synthesis_24h_10k_14_06_2022.RDS")

xgb_Synthesis_24h_predict_test <- predict(xgb_Synthesis_24h,Synthesis_24h_test)
xgb_Synthesis_24h_predict_test <- data.frame(xgb_Synthesis_24h_predict_test)

plot((xgb_Synthesis_24h_predict_test$xgb_Synthesis_24h_predict_test),(Synthesis_24h_test$log2N_per_min),asp = 1)

print(cor(log10(xgb_Synthesis_24h_predict_test$xgb_Synthesis_24h_predict_test),log10(Synthesis_24h_test$log2N_per_min),use = "pairwise.complete.obs",method = "pearson"))^2

gc()
rm(xgb_Synthesis_24h,Synthesis_24h_train,Synthesis_24h_test)













