library(plyr)
library(dplyr)
library(doMC)
library(randomForest)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)
library(ggpointdensity)
library(viridis)
library(xgboost)




###__________________________________________________________________________________________###
###____________________________________________Ex vivo_______________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")

CN_CD4_Tn_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_train",sep=";")
CN_CD4_Tn_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tn_test",sep=";")



## Training the xgb model ##
print("modeling CD4 TN CN")
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
                       nthread=36,
                       verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(xgb_CD4_Tn_CN)

saveRDS(xgb_CD4_Tn_CN, "~/Analysis/models/protein/xgb_CN_CD4_Tn_18_10_2021.RDS")
#xgb_CD4_Tn_CN <- readRDS("~/Analysis/models/protein/xgb_CN_CD4_Tn_18_10_2021.RDS")

xgb_CD4_Tn_CN_predict_test <- predict(xgb_CD4_Tn_CN,CN_CD4_Tn_test)
xgb_CD4_Tn_CN_predict_test <- data.frame(xgb_CD4_Tn_CN_predict_test)
print(lm(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test~CN_CD4_Tn_test$CN))

print(cor(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN, method = "pearson", use = "complete.obs")^2)
plot(xgb_CD4_Tn_CN_predict_test$xgb_CD4_Tn_CN_predict_test,CN_CD4_Tn_test$CN)




## Noise ## 

# 
# 
# CN_CD4_Tem_train$CN <- sample(CN_CD4_Tem_train$CN)
# 
# 
# 
# 
# ## Training the xgb model ##
# print("modeling CD4 Tem CN")
# registerDoMC(1)
# set.seed(12345)
# control <- trainControl(method="cv",
#                         #p=0.8,
#                         number=10,
#                         #repeats = 1,
#                         verboseIter = TRUE,
#                         allowParallel = TRUE,
#                         seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)
# 
# xgbGrid_CD4_Tem_CN_Noise <- expand.grid(nrounds = 1000,  
#                                   max_depth = 2,
#                                   colsample_bytree = 0.75,#seq(0.3, 0.5, length.out = 3),
#                                   eta = 0.05,
#                                   gamma=1,
#                                   min_child_weight = 0.9,
#                                   subsample = 0.75)
# 
# 
# start_time <- Sys.time()
# xgb_CD4_Tem_CN_Noise <- train(CN~.,
#                         data=CN_CD4_Tem_train,
#                         method="xgbTree",
#                         trControl=control,
#                         #metric="Rsquared",
#                         tuneGrid= xgbGrid_CD4_Tem_CN_Noise,
#                         na.action = na.omit,
#                         nthread=24,
#                         verbose = TRUE)
# 
# 
# end_time <- Sys.time()
# print(end_time - start_time) # Time difference of 
# 
# print(xgb_CD4_Tem_CN_Noise)
# 
# saveRDS(xgb_CD4_Tem_CN_Noise, "~/Analysis/models/protein/xgb_CN_Noise_CD4_Tem_22_10_2021.RDS")
xgb_CD4_Tem_CN_Noise <- readRDS("~/Analysis/models/protein/xgb_CN_Noise_CD4_Tem_22_10_2021.RDS")

xgb_CD4_Tem_CN_Noise_predict_test <- predict(xgb_CD4_Tem_CN_Noise,CN_CD4_Tem_test)
xgb_CD4_Tem_CN_Noise_predict_test <- data.frame(xgb_CD4_Tem_CN_Noise_predict_test)
print(lm(xgb_CD4_Tem_CN_Noise_predict_test$xgb_CD4_Tem_CN_Noise_predict_test~CN_CD4_Tem_test$CN))

print(cor(xgb_CD4_Tem_CN_Noise_predict_test$xgb_CD4_Tem_CN_Noise_predict_test,CN_CD4_Tem_test$CN, method = "pearson", use = "complete.obs")^2)

# 
# 
# 
