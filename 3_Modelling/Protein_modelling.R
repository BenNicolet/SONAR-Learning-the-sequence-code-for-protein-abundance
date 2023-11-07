# Serial mo

library(plyr)
library(dplyr)
library(doMC)
# library(biomaRt)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)
library(xgboost)
library(ggpointdensity)
library(viridis)
# library(readsxl)
# install.packages("biomaRt")

setwd("./")



###__________________________________________________________________________________________###
###___________________________________________lib_param______________________________________###
###__________________________________________________________________________________________###


## parameters ##
print("importing lib")
Sequence_parameters_protein <- read.delim("/DATA/users/b.nicolet/Cancer_cell_modeling/data/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")



###__________________________________________________________________________________________###
###____________________________________importing and prep-ing_________________________________###
###__________________________________________________________________________________________###


## prep data ##
# print("prep-ing data")
MS_immune_and_lines <- read.delim("./data/Cell_lines_and_immune_cells_CN_all_replicates_log10.csv", sep = ";", dec = ",") # importing the CN data 
# dim(MS_immune_and_lines)

# Here I make sure that the CN values are interpreted as numerical values 
MS_immune_and_lines[2:102] <- lapply(MS_immune_and_lines[2:102], function(x) {
  if(is.character(x)) as.numeric(x) else x
})


MS_immune_and_lines <- MS_immune_and_lines[MS_immune_and_lines$Entry %in% Sequence_parameters_protein$Entry,] # Here I select the CN rows that have a matching ID in the sequence feature library
# table(duplicated(MS_immune_and_lines$Entry)) # no dups

MS_immune_and_lines <- merge(MS_immune_and_lines,Sequence_parameters_protein,by.x="Entry",by.y="Entry", all.x=F) # merging data with SF library
rownames(MS_immune_and_lines) <- MS_immune_and_lines$Entry # adding rownames
MS_immune_and_lines$Entry <- NULL
# print(dim(MS_immune_and_lines)) # 10075  7225 # checking if size makes sense
# Here identify the last column containing CN data


## Feature selection ##
registerDoMC(6) # to allow parallel backend 

MS_immune_and_lines[is.na(MS_immune_and_lines)]=0
MS_immune_and_lines <- MS_immune_and_lines[, -nearZeroVar(MS_immune_and_lines, allowParallel = T, uniqueCut = 0.1)] # here I drop the feautres that have near-zero variance
print(dim(MS_immune_and_lines)) # 6013 non zero var features

gc()
registerDoMC(1) # Bring back to 1 core. This is very important!! XGBoost seems to paralellize things already. If core>1 then everything is run in double ( use nthread in caret instead to throttle the core # used)


## Modelling parameters ##
control <- caret::trainControl(method="cv", # Here I use a 5-fold cross-validation 
                               number=5,
                               verboseIter = TRUE, # talk to me
                               allowParallel = TRUE, # allowing things to run parallel
                               seeds = NULL) # NULL will set the seeds using a random set of integers

# xgbGrid allows you to find the optimal hyperparameters for your task. I already performed this and the optimal parameters are listed bellow. 
# For optimization, nrounds should be large, max depth >1, 0 >= colsample_bytree >= 1, eta range is typically 0.01-0.3, gamma is typically 0-10, min_child_weight is typically 0.4-1, subsample is typically 0.5-1
# Keep in mind that the hyperparameter optimization demand increases x2 with every new parameter. e.g. #cv * #nrounds * #max_depth * #colsample_bytree * #eta * #gamma * #min_child_weight * #subsample. 
# If you use a 5-fold CV and have 2 value

xgbGrid <- expand.grid(nrounds = 10000, # this will run the optimisation rounds 10k times. Rule of thumb nround > #features
                       max_depth = 6, # this is the number of levels allowed per tree. I tested many depth, but >6 for SONAR is good. The higer the value, the longer the modelling time.
                       colsample_bytree = 0.4,
                       eta = 0.05, # This controls the rate of prior knowledge for an iteration. typically 0.01-0.1
                       gamma=1, # This can be >1 too. 0 allows all SF to be considered and retained in the final model. gamma=1 prunes the nodes (and SFs) that contribute < gamma
                       min_child_weight = 0.9, 
                       subsample = 1) # here I want all of the data to be considered per iteration. lower values allow for faster modelling, but may come with costs to accuracy

set.seed(12345)# setting seed. If using Rmarkdown, set seed in EACH code chunks. 

## Serial modelling ##
cor_list <- list() # creating a holder list to store the correlation (test versus predicted)
Var_imp_list <- list() # creating a holder list to store the feature importance

{
  total_time <- Sys.time()
  for (i in 1:101) { # Here, change 101 for the last column which contains CN data.
    set.seed(12345) # just to be sure
    data <- cbind.data.frame(MS_immune_and_lines[i],MS_immune_and_lines[102:dim(MS_immune_and_lines)[2]]) # Here I take the i-th column and fuse it to the feature library. library starts 1 after the last CN data column.
    colnames(data)[1] <- "CN" # Giving a standard name
    
    data <- subset(data, data$CN>0) # getting rid of rows with 0s in CN data
    
    # split test/train sets
    train_row <- sample(1:nrow(data), 0.8 * nrow(data)) # Here I split a test and train set (80% recommended), to prevent data leakages. 
    data_train <- data[train_row,]
    data_test <- data[-train_row,]
    
    # Saving the test and train sets for future use.
    write.table(data_test,(paste0("./test_set_",colnames(MS_immune_and_lines[i]))),sep=";")
    write.table(data_train,(paste0("./train_set_",colnames(MS_immune_and_lines[i]))),sep=";")
    
    #
    
    ## Training the xgb model ##
    print(paste0("modeling ",colnames(MS_immune_and_lines[i]))) # Handy to report where the modelling is at. 
    start_model_time <- Sys.time()
    
    Model <- caret::train(CN~.,
                          data=data_train,
                          method="xgbTree",
                          trControl=control,# This points to the parameters define above
                          tuneGrid= xgbGrid,# This points to the parameters define above
                          na.action = na.omit, # This is pretty important, else the prediction of the test set may go funky.
                          nthread=48, # Set the "core" (actually thread) number, (total #core - 1)/2 works best (especially with intel virtualization with doubles the cores virutally)
                          verbose = TRUE) # talk to me 
    
    
    end_model_time <- Sys.time()
    print(paste0("Time difference of ", end_model_time - start_model_time)) # showing how long it took for this model to be made
    
    # Saving the model
    saveRDS(Model, (paste0("/DATA/users/b.nicolet/immune_cell_modelling/deeper_trees/models/XGB_model_",colnames(MS_immune_and_lines[i]),".RDS")))
    
    # predicting the held-out test set 
    predicted_test <- predict(Model,data_test)
    predicted_test <- data.frame("pred"=predicted_test)
    
    print(cor(predicted_test$pred,data_test$CN, method = "pearson", use = "complete.obs")^2) # This gives the Rsquared value
    cor_list[[colnames(MS_immune_and_lines[i])]] <- (cor(predicted_test$pred, data_test$CN, method = "pearson", use = "complete.obs")^2) # I aggregate the Rsquared values in this dataframe
    
    
    feature_imp <- data.frame(caret::varImp(Model,scale=F)$importance) # Extracting the feature importance of the models
    feature_imp$ID <- rownames(feature_imp)
    filename_var_imp <- paste0("/DATA/users/b.nicolet/immune_cell_modelling/deeper_trees/feature_imp/XGB_model_",colnames(MS_immune_and_lines[i]),".csv") 
    print(filename_var_imp)
    write.table(feature_imp, file = filename_var_imp, sep=";",row.names = F) # storing the feature importance in a separate file
    colnames(feature_imp)[1] <- colnames(MS_immune_and_lines[i])
    Var_imp_list[[colnames(MS_immune_and_lines[i])]] <- feature_imp # Aggregating the feature importance in one data frame for future export
    
    gc()# cleaning mem:
    
  }
  
  end_time <- Sys.time()
  end_time - total_time # How long did it take? 
}

# preparing the Rsquared performances for export
model_perf <- do.call(rbind,cor_list)
model_perf <- data.frame(model_perf) 

write.table(model_perf,file = "/DATA/users/b.nicolet/immune_cell_modelling/deeper_trees/model_perf_deep_trees_CN_immune_cell_replicates.csv",sep=";", row.names = T)

# preparing the feature importance df for export
feature_importance_all <- join_all(Var_imp_list, by = "ID")
write.table(feature_importance_all,file = "/DATA/users/b.nicolet/immune_cell_modelling/deeper_trees/feature_imp_deep_trees_CN_immune_cell_replicates.csv",sep=";", row.names = F)





