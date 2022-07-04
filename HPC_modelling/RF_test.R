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




###__________________________________________________________________________________________###
###____________________________________________Ex vivo_______________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")


CN_CD4_Tem_train <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_train",sep=";")
CN_CD4_Tem_test <- read.delim("~/Analysis/split_dataset/CN_CD4_Tem_test",sep=";")


dim(CN_CD4_Tem_train)
6128/3

## Training the rf model ##
print("modeling CD4 Tem CN")
registerDoMC(25)
set.seed(12345)
control <- trainControl(method="cv",
                        number=5,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

rfGrid_CD4_Tem_CN <- expand.grid(mtry=c(1500,2043,2500,3000,6128))

start_time <- Sys.time()
rf_CD4_Tem_CN <- train(CN~.,
                        data=CN_CD4_Tem_train,
                        method="rf",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= rfGrid_CD4_Tem_CN,
                        na.action = na.omit,
                        ntree=100,
                        verbose = TRUE)


end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CD4_Tem_CN)

saveRDS(rf_CD4_Tem_CN, "~/Analysis/models/protein/rf_CN_CD4_Tem_29_10_2021.RDS")
#rf_CD4_Tem_CN <- readRDS("~/Analysis/models/protein/rf_CN_CD4_Tem_18_10_2021.RDS")

rf_CD4_Tem_CN_predict_test <- predict(rf_CD4_Tem_CN,CN_CD4_Tem_test)
rf_CD4_Tem_CN_predict_test <- data.frame(rf_CD4_Tem_CN_predict_test)
print(lm(rf_CD4_Tem_CN_predict_test$rf_CD4_Tem_CN_predict_test~CN_CD4_Tem_test$CN))+ 



gc()
rm(rf_CD4_Tem_CN,CN_CD4_Tem)





