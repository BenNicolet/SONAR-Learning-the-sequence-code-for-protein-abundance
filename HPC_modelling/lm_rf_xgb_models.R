library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)


# 
# 
# ###__________________________________________________________________________________________###
# ###____________________________________________Ex vivo_______________________________________###
# ###__________________________________________________________________________________________###
# 
# setwd("~/Analysis/")
# 
# 
# ###__________________________________________________________________________________________###
# ###___________________________________________CD8_Tn_______________________________________###
# ###__________________________________________________________________________________________###
# 
# ## prep data ##
# print("preping data for CD8 Tn")
# TPM_CD8_Tn_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_train",sep=";")
# TPM_CD8_Tn_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_test",sep=";")
# 
# 
# 
# 
# gc()
# 
# 
# ## Training the lm model for Tn 0h ##
# registerDoMC(20)
# set.seed(12345)
# control <- trainControl(method="repeatedcv",
#                         number=10,
#                         repeats = 2,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# # lmGrid_CD8Tn <- expand.grid(nrounds = 10000,  # this is n_estimators in the python code above
# #                              max_depth = 6,
# #                              colsample_bytree = 0.4,#seq(0.3, 0.5, length.out = 3),
# #                              ## The values below are default values in the sklearn-api.
# #                              eta = 0.05,
# #                              gamma=0,
# #                              min_child_weight = 0.9,
# #                              subsample = 1)
# 
# start_time <- Sys.time()
# lm_TPM_CD8_Tn <- train(TPM~.,
#                         data=TPM_CD8_Tn_train,
#                         method="lm",
#                         trControl=control,
#                         #metric="Rsquared",
#                         #tuneGrid= lmGrid_CD8Tn,
#                         na.action = na.omit,
#                         #nthread=94,
#                         verbose = TRUE)
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of
# print(lm_TPM_CD8_Tn)
# 
# saveRDS(lm_TPM_CD8_Tn, "~/Analysis/models/lm_TPM_CD8_Tn_08_11_2021.RDS")
# #lm_TPM_CD8_Tn <- readRDS("~/Analysis/models/lm_TPM_CD8_Tn_08_11_2021.RDS")
# 
lm_TPM_CD8_Tn_predict_test <- predict(lm_TPM_CD8_Tn,TPM_CD8_Tn_test)
lm_TPM_CD8_Tn_predict_test <- data.frame(lm_TPM_CD8_Tn_predict_test)
print(lm(lm_TPM_CD8_Tn_predict_test$lm_TPM_CD8_Tn_predict_test~TPM_CD8_Tn_test$TPM))
plot(lm_TPM_CD8_Tn_predict_test$lm_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM,xlim = c(-2,4),ylim = c(-2,4))+
  abline(lm(lm_TPM_CD8_Tn_predict_test$lm_TPM_CD8_Tn_predict_test~TPM_CD8_Tn_test$TPM),col='red')


cor(lm_TPM_CD8_Tn_predict_test$lm_TPM_CD8_Tn_predict_test, TPM_CD8_Tn_test$TPM,  method = "pearson", use = "complete.obs")^2

# gc()
# rm(lm_TPM_CD8_Tn,TPM_CD8_Tn_train,TPM_CD8_Tn_test)
# 
# 




###__________________________________________________________________________________________###
###____________________________________________Ex vivo_______________________________________###
###__________________________________________________________________________________________###

setwd("~/Analysis/")


###__________________________________________________________________________________________###
###___________________________________________CD8_Tn_______________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data for CD8 Tn")
TPM_CD8_Tn_train <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_train",sep=";")
TPM_CD8_Tn_test <- read.delim("~/Analysis/split_dataset/TPM_CD8_Tn_test",sep=";")




gc()


## Training the rf model for Tn 0h ##
registerDoMC(20)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

rfGrid_CD8Tn <- expand.grid(mtry = round(6128/3))

start_time <- Sys.time()
rf_TPM_CD8_Tn <- train(TPM~.,
                       data=TPM_CD8_Tn_train,
                       method="rf",
                       trControl=control,
                       ntree=1000,
                       #metric="Rsquared",
                       tuneGrid= rfGrid_CD8Tn,
                       na.action = na.omit,
                       #nthread=94,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of
print(rf_TPM_CD8_Tn)

saveRDS(rf_TPM_CD8_Tn, "~/Analysis/models/rf_TPM_CD8_Tn_08_11_2021.RDS")
rf_TPM_CD8_Tn <- readRDS("~/Analysis/models/rf_TPM_CD8_Tn_08_11_2021.RDS")

rf_TPM_CD8_Tn_predict_test <- predict(rf_TPM_CD8_Tn,TPM_CD8_Tn_test)
rf_TPM_CD8_Tn_predict_test <- data.frame(rf_TPM_CD8_Tn_predict_test)
print(rf(rf_TPM_CD8_Tn_predict_test$rf_TPM_CD8_Tn_predict_test~TPM_CD8_Tn_test$TPM))
plot(rf_TPM_CD8_Tn_predict_test$rf_TPM_CD8_Tn_predict_test,TPM_CD8_Tn_test$TPM)

cor(rf_TPM_CD8_Tn_predict_test$rf_TPM_CD8_Tn_predict_test, TPM_CD8_Tn_test$TPM,  method = "pearson", use = "complete.obs")^2


gc()
rm(rf_TPM_CD8_Tn,TPM_CD8_Tn_train,TPM_CD8_Tn_test)



