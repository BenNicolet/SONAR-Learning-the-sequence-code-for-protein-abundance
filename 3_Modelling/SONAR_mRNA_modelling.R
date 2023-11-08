## SONAR mRNA modelling script ##
## Ben Nicolet ##


library(plyr)
library(dplyr)
library(doMC)
library(tidyverse)
library(caret)
library(e1071)
library(xgboost)
library(doMC)


setwd("./")


###__________________________________________________________________________________________###
###__________________________________________SF library______________________________________###
###__________________________________________________________________________________________###

## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("./RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB.csv", sep=";", dec=",")
Sequence_parameters_RNA$gene_name <- NULL
Sequence_parameters_RNA$external_gene_name <- NULL


tx2gene <- data.frame("tx.id"=Sequence_parameters_RNA$tx.id,
                      "ensembl_gene_id"=Sequence_parameters_RNA$ensembl_gene_id)
Sequence_parameters_RNA <- as.data.frame(Sequence_parameters_RNA)
Sequence_parameters_RNA[3:dim(Sequence_parameters_RNA)[2]] <- lapply(Sequence_parameters_RNA[3:dim(Sequence_parameters_RNA)[2]], function(x) {
  if(is.character(x)) as.numeric(x) else x
})



###__________________________________________________________________________________________###
###____________________________________importing and preping_________________________________###
###__________________________________________________________________________________________###



## prep data ##
print("prep-ing data")
TPM_immune_and_lines <- read.delim("./data/TPM_immune_and_cell_lines_TPM-2.csv", sep = ";", dec = ",") # getting TPMs and dropping a few columns
TPM_immune_and_lines$external_gene_name <- NULL
rownames(TPM_immune_and_lines) <- TPM_immune_and_lines$ID
TPM_immune_and_lines$ID <- NULL

## Modelling parameters ## 
# See SONAR protein script for details ##
control <- caret::trainControl(method="cv",
                               number=5,
                               verboseIter = TRUE,
                               allowParallel = TRUE,
                               seeds = NULL) 

xgbGrid <- expand.grid(nrounds = 10000, 
                       max_depth = 6,
                       colsample_bytree = 0.4,
                       eta = 0.05,
                       gamma=1,
                       min_child_weight = 0.9,
                       subsample = 1)

set.seed(12345)



## Preparation for and modelling ##
{
  total_time <- Sys.time()
  for (i in 1:88) { # here should be 1:#_of_samples
    
    set.seed(12345)
    start_model_time <- Sys.time()
    # Here I split 1 column from the data (i-th column) and merge it with the feature library
    data <- data.frame("ID"=rownames(TPM_immune_and_lines), TPM_immune_and_lines[i])
    colnames(data)[2] <- "TPM"
    data <- subset(data, data$TPM>=0.1)
    data <- merge(data, Sequence_parameters_RNA, by.x="ID",by.y="tx.id", all.x=F)



    print( paste0("Getting SFs per genes for: ",colnames(TPM_immune_and_lines)[i]))
    # As we have TPM per transcript, I average the corresponding feautre info 
    SF_lib <- data %>% group_by(ensembl_gene_id) %>% summarise(across(.fns =  mean, na.rm = TRUE))
    SF_lib$TPM <- NULL
    print(paste0("SF_lib dims: ",dim(SF_lib)))



    # Here I drop transcripts with less than 0.1TPM and sum up TPM per gene. Of note, this step can be dropped if there is interest in tanscript level modelling (be aware of the high computing demand in that case...)
    print( paste0("Getting TPM per genes for: ",colnames(TPM_immune_and_lines)[i]))
    tx_list <- data.frame("ensembl_gene_id"=data$ensembl_gene_id, "TPM"=data$TPM)
    tx_list <- subset(tx_list, tx_list$TPM>=0.1)
    tx_list <- tx_list %>% group_by(ensembl_gene_id) %>% summarise(across(.fns =  sum, na.rm = TRUE))
    print(paste0("tx_list dims: ",dim(tx_list)))
    tx_list$TPM <- log10(tx_list$TPM)


    # Merging data + feature (both per genes) together
    data <- merge(tx_list, SF_lib, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
    rownames(data) <- data$ensembl_gene_id
    data$ensembl_gene_id <- NULL
    data$tx.id <- NULL
    print(paste0("agregated data are ", dim(data)[1]," by ",dim(data)[2]))

    # Cleaning up the features with zero and near-zero variance 
    data[is.na(data)]=0
    data <- data[, -nearZeroVar(data, allowParallel = T, uniqueCut = 0.1)]
    print(paste0("Non-zero var agregated data are ", dim(data)[1]," by ",dim(data)[2]))

    gc()
    registerDoMC(1)



    # split test/train sets
    train_row <- sample(1:nrow(data), 0.8 * nrow(data))
    data_train <- data[train_row,]
    data_test <- data[-train_row,]
    
    # Saving split sets for future use
    write.table(data_test,(paste0("./split_datasets/test_set_",colnames(TPM_immune_and_lines[i]))),sep=";")
    write.table(data_train,(paste0("./split_datasets/train_set_",colnames(TPM_immune_and_lines[i]))),sep=";")
    
    # If test/train sets have already been defined: 
    # data_test <- read.delim((paste0("./split_datasets/test_set_",colnames(TPM_immune_and_lines[i]))),sep=";")
    # data_train <- read.delim((paste0("./split_datasets/train_set_",colnames(TPM_immune_and_lines[i]))),sep=";")
    
    
    
    
    
    ## Training the  model ##
    registerDoMC(1) # Do no increase this. Use nthread option of caret!!!!
    print(paste0("modeling ",colnames(TPM_immune_and_lines[i])))
    
    # Training 
    Model <- caret::train(TPM~.,
                          data=data_train,
                          method="xgbTree",
                          trControl=control,
                          tuneGrid= xgbGrid,
                          na.action = na.omit,
                          nthread=48,
                          verbose = TRUE)
    
    
    end_model_time <- Sys.time()
    print(paste0("Time difference of ", end_model_time - start_model_time))
    # Saving 
    saveRDS(Model, (paste0("./XGB_model_",colnames(TPM_immune_and_lines[i]),".RDS")))
    
    # predicting the held out testset
    predicted_test <- predict(Model,data_test)
    predicted_test <- data.frame("pred"=predicted_test)
    # and comparing it to the ground truth to get an Rsquared
    print(paste0("Rsquare = ",cor(predicted_test$pred,data_test$TPM, method = "pearson", use = "complete.obs")^2))
    print(paste0("Modelling at : ",((i/88)*100),"%"))
    gc()
    
  }
  
  end_time <- Sys.time()
  end_time - total_time
}








#### Feature importance, model metrics ####


## importing TPM to get colnames ##
TPM_immune_and_lines <- read.delim("./data/TPM_immune_and_cell_lines_TPM-2.csv", sep = ";", dec = ",")
TPM_immune_and_lines$external_gene_name <- NULL
rownames(TPM_immune_and_lines) <- TPM_immune_and_lines$ID
TPM_immune_and_lines$ID <- NULL


#___________________________________#
### Feature importance extraction ###
#___________________________________#


# 4) make a for loop to: 
# - import model
# - import test
# - extract correlation
# - extract feature importance


cor_list <- list()
Var_imp_list <- list()

{
  total_time <- Sys.time()
  for (i in 1:88) {
    
    # Import model
    model_filename <- paste0("./models_gamma=1/XGB_model_",colnames(TPM_immune_and_lines[i]),".RDS")
    print(model_filename)
    model <- readRDS(model_filename)
    
    
    
    ## feature imp ##
    feature_imp <- data.frame(caret::varImp(model,scale=F)$importance)
    feature_imp$ID <- rownames(feature_imp)
    filename_var_imp <- paste0("./feature_imp_gamma1/XGB_model_",colnames(TPM_immune_and_lines[i]),".csv")
    print(filename_var_imp)
    write.table(feature_imp, file = filename_var_imp, sep=";",row.names = F)
    colnames(feature_imp)[1] <- colnames(TPM_immune_and_lines[i])
    Var_imp_list[[colnames(TPM_immune_and_lines[i])]] <- feature_imp
    
    
    
    ## importing test sets ##
    filename_test <- paste0("./split_datasets/test_set_",colnames(TPM_immune_and_lines[i]))
    print(filename_test)
    data_test <- read.delim(filename_test,sep=";")
    
    
    
    ## Rsq ##
    predicted_test <- predict(model,data_test)
    predicted_test <- data.frame("pred"=predicted_test)
    
    print(cor(predicted_test$pred,data_test$TPM, method = "pearson", use = "complete.obs")^2)
    cor_list[[colnames(TPM_immune_and_lines[i])]] <- (cor(predicted_test$pred,data_test$TPM, method = "pearson", use = "complete.obs")^2)
    
    # cleaning mem:
    gc()
    
  }
  
  end_time <- Sys.time()
  print(end_time - total_time)
  
}


model_perf <- do.call(rbind,cor_list)
model_perf <- data.frame(model_perf)
feature_importance_all <- join_all(Var_imp_list, by = "ID")


write.table(model_perf,file = "./model_performance_gamma1_RNA_immune_cell_replicates.csv",sep=";", row.names = T)
write.table(feature_importance_all,file = "./feature_importance_gamma1_RNA_immune_cell_replicates.csv",sep=";", row.names = F)








