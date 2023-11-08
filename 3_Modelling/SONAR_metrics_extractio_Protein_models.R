## Formatting script for SONAR output ##
## Ben Nicolet ##
## 07-10-2022 ##

library(plyr)

# This script should be run if the main modelling script bugs halfway, as it will export all metrics needed for analysis, minus the modelling task

# make a for loop to: 
# - import model
# - import test
# - extract correlation
# - extract feature importance

setwd("./")

## prep data ##
print("prep-ing data")
MS_immune_and_lines <- read.delim("./Cell_lines_and_immune_cells_CN_all_replicates_log10.csv", sep = ";", dec = ",")
# dim(MS_immune_and_lines)

MS_immune_and_lines[2:102] <- lapply(MS_immune_and_lines[2:102], function(x) {
  if(is.character(x)) as.numeric(x) else x
})



rownames(MS_immune_and_lines) <- MS_immune_and_lines$Entry
MS_immune_and_lines$Entry <- NULL
print(dim(MS_immune_and_lines)) # 10075  7225
# data of cells stop at 101





#___________________#
### FI extraction ###
#___________________#


# 4) make a for loop to: 
# - import model
# - import test
# - extract correlation
# - extract feature importance


cor_list <- list()
Var_imp_list <- list()

{
  total_time <- Sys.time()
  for (i in 1:101) {
    
    # Import model
    model_filename <- paste0("./models/XGB_model_",colnames(MS_immune_and_lines[i]),".RDS")
    print(model_filename)
    model <- readRDS(model_filename)
    
    
    
    ## feature imp ##
    feature_imp <- data.frame(caret::varImp(model,scale=F)$importance)
    feature_imp$ID <- rownames(feature_imp)
    filename_var_imp <- paste0("./feature_imp/XGB_model_",colnames(MS_immune_and_lines[i]),".csv")
    print(filename_var_imp)
    write.table(feature_imp, file = filename_var_imp, sep=";",row.names = F)
    colnames(feature_imp)[1] <- colnames(MS_immune_and_lines[i])
    Var_imp_list[[colnames(MS_immune_and_lines[i])]] <- feature_imp
    
    
    
    ## importing test sets ##
    filename_test <- paste0("./split_datasets/test_set_",colnames(MS_immune_and_lines[i]))
    print(filename_test)
    data_test <- read.delim(filename_test,sep=";")
    
    
    
    ## Rsq ##
    predicted_test <- predict(model,data_test)
    predicted_test <- data.frame("pred"=predicted_test)
    
    print(cor(predicted_test$pred,data_test$CN, method = "pearson", use = "complete.obs")^2)
    cor_list[[colnames(MS_immune_and_lines[i])]] <- (cor(predicted_test$pred,data_test$CN, method = "pearson", use = "complete.obs")^2)
    
    # cleaning mem:
    gc()
    
  }
  
  end_time <- Sys.time()
  print(end_time - total_time)
  
}


model_perf <- do.call(rbind,cor_list)
model_perf <- data.frame(model_perf)
feature_importance_all <- join_all(Var_imp_list, by = "ID")


write.table(model_perf,file = "./model_performance_CN_immune_cell_replicates.csv",sep=";", row.names = T)
write.table(feature_importance_all,file = "./feature_importance_CN_immune_cell_replicates.csv",sep=";", row.names = F)


