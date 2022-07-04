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


## parameters ##
Sequence_parameters_protein <- read.delim("~/Analysis/libraries/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv", sep="\t", dec=".")





###__________________________________________________________________________________________###
###____________________________________________CD8_Tn________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Tn <- read.delim("~/Analysis/data/protein/CD8_Tn_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tn <- subset(CN_CD8_Tn,CN_CD8_Tn$CN>0)
dim(CN_CD8_Tn)

CN_CD8_Tn <- merge(CN_CD8_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tn[is.na(CN_CD8_Tn)]= 0
print(dim(CN_CD8_Tn)) # 7984 7126


rownames(CN_CD8_Tn) <- CN_CD8_Tn$ID
CN_CD8_Tn$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tn[is.na(CN_CD8_Tn)]=0
CN_CD8_Tn <- CN_CD8_Tn[, -nearZeroVar(CN_CD8_Tn, allowParallel = T, uniqueCut = dim(CN_CD8_Tn)[1]*0.01)]
print(dim(CN_CD8_Tn)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD8_Tn <- train(CN~.,
                      data=CN_CD8_Tn,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Tn)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD8_Tn)
print(rf_CN_CD8_Tn$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 919
# Mean of squared residuals: 0.7148345
# % Var explained: 39.57

saveRDS(rf_CN_CD8_Tn, "~/Analysis/models/protein/rf_CN_CD8_Tn_09_09_2021.RDS")


gc()
rm(rf_CN_CD8_Tn,CN_CD8_Tn)


###__________________________________________________________________________________________###
###___________________________________________CD8_Tcm________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Tcm <- read.delim("~/Analysis/data/protein/CD8_Tcm_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tcm <- subset(CN_CD8_Tcm,CN_CD8_Tcm$CN>0)

CN_CD8_Tcm <- merge(CN_CD8_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tcm[is.na(CN_CD8_Tcm)]= 0
print(dim(CN_CD8_Tcm)) # 7984 7126


rownames(CN_CD8_Tcm) <- CN_CD8_Tcm$ID
CN_CD8_Tcm$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tcm[is.na(CN_CD8_Tcm)]=0
CN_CD8_Tcm <- CN_CD8_Tcm[, -nearZeroVar(CN_CD8_Tcm, allowParallel = T, uniqueCut = dim(CN_CD8_Tcm)[1]*0.01)]
print(dim(CN_CD8_Tcm)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD8_Tcm <- train(CN~.,
                      data=CN_CD8_Tcm,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Tcm)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD8_Tcm)
print(rf_CN_CD8_Tcm$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 916
# Mean of squared residuals: 0.6885122
# % Var explained: 39.39

saveRDS(rf_CN_CD8_Tcm, "~/Analysis/models/protein/rf_CN_CD8_Tcm_09_09_2021.RDS")


gc()
rm(rf_CN_CD8_Tcm,CN_CD8_Tcm)





###__________________________________________________________________________________________###
###___________________________________________CD8_Tem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Tem <- read.delim("~/Analysis/data/protein/CD8_Tem_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Tem <- subset(CN_CD8_Tem,CN_CD8_Tem$CN>0)

CN_CD8_Tem <- merge(CN_CD8_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Tem[is.na(CN_CD8_Tem)]= 0
print(dim(CN_CD8_Tem)) # 7984 7126


rownames(CN_CD8_Tem) <- CN_CD8_Tem$ID
CN_CD8_Tem$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Tem[is.na(CN_CD8_Tem)]=0
CN_CD8_Tem <- CN_CD8_Tem[, -nearZeroVar(CN_CD8_Tem, allowParallel = T, uniqueCut = dim(CN_CD8_Tem)[1]*0.01)]
print(dim(CN_CD8_Tem)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD8_Tem <- train(CN~.,
                       data=CN_CD8_Tem,
                       method="rf",
                       ntree=5000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Tem)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD8_Tem)
print(rf_CN_CD8_Tem$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 915
# Mean of squared residuals: 0.6888405
# % Var explained: 39.27

saveRDS(rf_CN_CD8_Tem, "~/Analysis/models/protein/rf_CN_CD8_Tem_09_09_2021.RDS")

gc()
rm(rf_CN_CD8_Tem,CN_CD8_Tem)





###__________________________________________________________________________________________###
###___________________________________________CD8_Teff________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD8_Teff <- read.delim("~/Analysis/data/protein/CD8_Teff_CN_log10.csv",sep=";",dec = ",")
CN_CD8_Teff <- subset(CN_CD8_Teff,CN_CD8_Teff$CN>0)

CN_CD8_Teff <- merge(CN_CD8_Teff,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD8_Teff[is.na(CN_CD8_Teff)]= 0
print(dim(CN_CD8_Teff)) # 7984 7126


rownames(CN_CD8_Teff) <- CN_CD8_Teff$ID
CN_CD8_Teff$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD8_Teff[is.na(CN_CD8_Teff)]=0
CN_CD8_Teff <- CN_CD8_Teff[, -nearZeroVar(CN_CD8_Teff, allowParallel = T, uniqueCut = dim(CN_CD8_Teff)[1]*0.01)]
print(dim(CN_CD8_Teff)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD8_Teff <- train(CN~.,
                       data=CN_CD8_Teff,
                       method="rf",
                       ntree=5000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(CN_CD8_Teff)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD8_Teff)
print(rf_CN_CD8_Teff$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 915
# Mean of squared residuals: 0.7071294
# % Var explained: 36.6

saveRDS(rf_CN_CD8_Teff, "~/Analysis/models/protein/rf_CN_CD8_Teff_09_09_2021.RDS")

gc()
rm(rf_CN_CD8_Teff,CN_CD8_Teff)










###__________________________________________________________________________________________###
###____________________________________________CD4_Tn________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tn <- read.delim("~/Analysis/data/protein/CD4_Tn_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tn <- subset(CN_CD4_Tn,CN_CD4_Tn$CN>0)
dim(CN_CD4_Tn)

CN_CD4_Tn <- merge(CN_CD4_Tn,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tn[is.na(CN_CD4_Tn)]= 0
print(dim(CN_CD4_Tn)) # 7984 7126


rownames(CN_CD4_Tn) <- CN_CD4_Tn$ID
CN_CD4_Tn$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tn[is.na(CN_CD4_Tn)]=0
CN_CD4_Tn <- CN_CD4_Tn[, -nearZeroVar(CN_CD4_Tn, allowParallel = T, uniqueCut = dim(CN_CD4_Tn)[1]*0.01)]
print(dim(CN_CD4_Tn)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD4_Tn <- train(CN~.,
                      data=CN_CD4_Tn,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_CD4_Tn)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD4_Tn)
print(rf_CN_CD4_Tn$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 916
# Mean of squared residuals: 0.6751644
# % Var explained: 38.99

saveRDS(rf_CN_CD4_Tn, "~/Analysis/models/protein/rf_CN_CD4_Tn_09_09_2021.RDS")


gc()
rm(rf_CN_CD4_Tn,CN_CD4_Tn)


###__________________________________________________________________________________________###
###___________________________________________CD4_Tcm________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tcm <- read.delim("~/Analysis/data/protein/CD4_Tcm_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tcm <- subset(CN_CD4_Tcm,CN_CD4_Tcm$CN>0)

CN_CD4_Tcm <- merge(CN_CD4_Tcm,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tcm[is.na(CN_CD4_Tcm)]= 0
print(dim(CN_CD4_Tcm)) # 7984 7126


rownames(CN_CD4_Tcm) <- CN_CD4_Tcm$ID
CN_CD4_Tcm$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tcm[is.na(CN_CD4_Tcm)]=0
CN_CD4_Tcm <- CN_CD4_Tcm[, -nearZeroVar(CN_CD4_Tcm, allowParallel = T, uniqueCut = dim(CN_CD4_Tcm)[1]*0.01)]
print(dim(CN_CD4_Tcm)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD4_Tcm <- train(CN~.,
                       data=CN_CD4_Tcm,
                       method="rf",
                       ntree=5000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(CN_CD4_Tcm)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD4_Tcm)
print(rf_CN_CD4_Tcm$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 917
# Mean of squared residuals: 0.6932329
# % Var explained: 39.24

saveRDS(rf_CN_CD4_Tcm, "~/Analysis/models/protein/rf_CN_CD4_Tcm_09_09_2021.RDS")


gc()
rm(rf_CN_CD4_Tcm,CN_CD4_Tcm)





###__________________________________________________________________________________________###
###___________________________________________CD4_Tem________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Tem <- read.delim("~/Analysis/data/protein/CD4_Tem_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Tem <- subset(CN_CD4_Tem,CN_CD4_Tem$CN>0)

CN_CD4_Tem <- merge(CN_CD4_Tem,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Tem[is.na(CN_CD4_Tem)]= 0
print(dim(CN_CD4_Tem)) # 7984 7126


rownames(CN_CD4_Tem) <- CN_CD4_Tem$ID
CN_CD4_Tem$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Tem[is.na(CN_CD4_Tem)]=0
CN_CD4_Tem <- CN_CD4_Tem[, -nearZeroVar(CN_CD4_Tem, allowParallel = T, uniqueCut = dim(CN_CD4_Tem)[1]*0.01)]
print(dim(CN_CD4_Tem)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD4_Tem <- train(CN~.,
                       data=CN_CD4_Tem,
                       method="rf",
                       ntree=5000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(CN_CD4_Tem)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD4_Tem)
print(rf_CN_CD4_Tem$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 918
# Mean of squared residuals: 0.6876779
# % Var explained: 39.12

saveRDS(rf_CN_CD4_Tem, "~/Analysis/models/protein/rf_CN_CD4_Tem_09_09_2021.RDS")

gc()
rm(rf_CN_CD4_Tem,CN_CD4_Tem)





###__________________________________________________________________________________________###
###___________________________________________CD4_Teff________________________________________###
###__________________________________________________________________________________________###


## prep data ##
CN_CD4_Teff <- read.delim("~/Analysis/data/protein/CD4_Teff_CN_log10.csv",sep=";",dec = ",")
CN_CD4_Teff <- subset(CN_CD4_Teff,CN_CD4_Teff$CN>0)

CN_CD4_Teff <- merge(CN_CD4_Teff,Sequence_parameters_protein,by.x="ID",by.y="Entry", all.x=F)
CN_CD4_Teff[is.na(CN_CD4_Teff)]= 0
print(dim(CN_CD4_Teff)) # 7984 7126


rownames(CN_CD4_Teff) <- CN_CD4_Teff$ID
CN_CD4_Teff$ID <- NULL


## Feature selection ##
registerDoMC(6)

CN_CD4_Teff[is.na(CN_CD4_Teff)]=0
CN_CD4_Teff <- CN_CD4_Teff[, -nearZeroVar(CN_CD4_Teff, allowParallel = T, uniqueCut = dim(CN_CD4_Teff)[1]*0.01)]
print(dim(CN_CD4_Teff)) # 7984 2759

gc()


## Training the RF model ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_CD4_Teff <- train(CN~.,
                        data=CN_CD4_Teff,
                        method="rf",
                        ntree=5000,
                        trControl=control,
                        metric="Rsquared",
                        tuneGrid= data.frame(mtry = round((ncol(CN_CD4_Teff)-1)/3)),
                        na.action = na.omit,
                        importance = TRUE,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 

print(rf_CN_CD4_Teff)
print(rf_CN_CD4_Teff$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 918
# Mean of squared residuals: 0.6740351
# % Var explained: 38.15

saveRDS(rf_CN_CD4_Teff, "~/Analysis/models/protein/rf_CN_CD4_Teff_09_09_2021.RDS")

gc()
rm(rf_CN_CD4_Teff,CN_CD4_Teff)

















###__________________________________________________________________________________________###
###____________________________________T cell activation_____________________________________###
###__________________________________________________________________________________________###

ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# attributes_ens <- listAttributes(ensembl)
# View(attributes_ens)
IDs_genenames_coding <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","uniprotswissprot","gene_biotype"), mart = ensembl, filters = "with_ccds", values = TRUE)
#IDs_genenames_coding <- IDs_genenames_coding[!duplicated(IDs_genenames_coding$uniprotswissprot),]


## parameters ##
print("importing library")
Sequence_parameters_Geiger <- read.delim("libraries/Sequence_parameters_for_Geiger_data_prot_libv2_09092021.csv",sep = ";",dec=",")




###__________________________________________________________________________________________###
###___________________________________________CD4_T0h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
print("importing and preparing data")
CN_0h <- read.delim("data/protein/CD4_Tn_0h_CN_log10.csv",sep = ";",dec = ",")
CN_0h <- merge(CN_0h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(CN_0h)) # 6993 6803
rownames(CN_0h) <- CN_0h$ID
CN_0h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_0h[is.na(CN_0h)]=0
CN_0h <- CN_0h[, -nearZeroVar(CN_0h, allowParallel = T, uniqueCut = dim(CN_0h)[1]*0.01)]
print(dim(CN_0h)) # 13729  2682

gc()


## Training the RF model for Tn 0h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_0h <- train(CN~.,
                         data=CN_0h,
                         method="rf",
                         ntree=5000,
                         trControl=control,
                         metric="Rsquared",
                         tuneGrid= data.frame(mtry = round((ncol(CN_0h)-1)/3)),
                         na.action = na.omit,
                         importance = TRUE,
                         verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(rf_CN_Tn_0h)
print(rf_CN_Tn_0h$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 929
# Mean of squared residuals: 0.6024929
# % Var explained: 42.73

saveRDS(rf_CN_Tn_0h, "models/protein/rf_CN_Tn_0h_Lib_V2_11092021.RDS")

gc()
rm(rf_CN_Tn_0h,CN_0h)



###__________________________________________________________________________________________###
###___________________________________________CD4_T6h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
print("importing and preparing data")
CN_6h <- read.delim("data/protein/CD4_Tn_6h_CN_log10.csv",sep = ";",dec = ",")
CN_6h <- merge(CN_6h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(CN_6h)) # 6993 6803
rownames(CN_6h) <- CN_6h$ID
CN_6h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_6h[is.na(CN_6h)]=0
CN_6h <- CN_6h[, -nearZeroVar(CN_6h, allowParallel = T, uniqueCut = dim(CN_6h)[1]*0.01)]
print(dim(CN_6h)) # 13729  2682

gc()


## Training the RF model for Tn 6h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_6h <- train(CN~.,
                     data=CN_6h,
                     method="rf",
                     ntree=5000,
                     trControl=control,
                     metric="Rsquared",
                     tuneGrid= data.frame(mtry = round((ncol(CN_6h)-1)/3)),
                     na.action = na.omit,
                     importance = TRUE,
                     verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(rf_CN_Tn_6h)
print(rf_CN_Tn_6h$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 931
# Mean of squared residuals: 0.5045108
# % Var explained: 45.97

saveRDS(rf_CN_Tn_6h, "models/protein/rf_CN_Tn_6h_Lib_V2_11092021.RDS")

gc()
rm(rf_CN_Tn_6h,CN_6h)







###__________________________________________________________________________________________###
###___________________________________________CD4_T12h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
print("importing and preparing data")
CN_12h <- read.delim("data/protein/CD4_Tn_12h_CN_log10.csv",sep = ";",dec = ",")
CN_12h <- merge(CN_12h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(CN_12h)) # 6993 6803
rownames(CN_12h) <- CN_12h$ID
CN_12h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_12h[is.na(CN_12h)]=0
CN_12h <- CN_12h[, -nearZeroVar(CN_12h, allowParallel = T, uniqueCut = dim(CN_12h)[1]*0.01)]
print(dim(CN_12h)) # 13729  2682

gc()


## Training the RF model for Tn 12h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_12h <- train(CN~.,
                     data=CN_12h,
                     method="rf",
                     ntree=5000,
                     trControl=control,
                     metric="Rsquared",
                     tuneGrid= data.frame(mtry = round((ncol(CN_12h)-1)/3)),
                     na.action = na.omit,
                     importance = TRUE,
                     verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(rf_CN_Tn_12h)
print(rf_CN_Tn_12h$finalModel)

# Number of trees: 5000
# No. of variables tried at each split: 931
# Mean of squared residuals: 0.4798387
# % Var explained: 47.81

saveRDS(rf_CN_Tn_12h, "models/protein/rf_CN_Tn_12h_Lib_V2_11092021.RDS")

gc()
rm(rf_CN_Tn_12h,CN_12h)








###__________________________________________________________________________________________###
###___________________________________________CD4_T24h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
print("importing and preparing data")
CN_24h <- read.delim("data/protein/CD4_Tn_24h_CN_log10.csv",sep = ";",dec = ",")
CN_24h <- merge(CN_24h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(CN_24h)) # 6993 6803
rownames(CN_24h) <- CN_24h$ID
CN_24h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_24h[is.na(CN_24h)]=0
CN_24h <- CN_24h[, -nearZeroVar(CN_24h, allowParallel = T, uniqueCut = dim(CN_24h)[1]*0.01)]
print(dim(CN_24h)) # 13729  2682

gc()


## Training the RF model for Tn 24h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_24h <- train(CN~.,
                      data=CN_24h,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_24h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 
print(rf_CN_Tn_24h)
print(rf_CN_Tn_24h$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 930
# Mean of squared residuals: 0.5023284
# % Var explained: 47.45

saveRDS(rf_CN_Tn_24h, "models/protein/rf_CN_Tn_24h_Lib_V2_11092021.RDS")

gc()
rm(rf_CN_Tn_24h,CN_24h)






###__________________________________________________________________________________________###
###___________________________________________CD4_T48h________________________________________###
###__________________________________________________________________________________________###


## prep data ##
print("importing and preparing data")
CN_48h <- read.delim("data/protein/CD4_Tn_48h_CN_log10.csv",sep = ";",dec = ",")
CN_48h <- merge(CN_48h,Sequence_parameters_Geiger,by="ID",all.x=F)
print(dim(CN_48h)) # 6993 6803
rownames(CN_48h) <- CN_48h$ID
CN_48h$ID <- NULL


## Feature selection ##
print("Feature selection")
registerDoMC(4)

CN_48h[is.na(CN_48h)]=0
CN_48h <- CN_48h[, -nearZeroVar(CN_48h, allowParallel = T, uniqueCut = dim(CN_48h)[1]*0.01)]
print(dim(CN_48h)) # 13729  2682

gc()


## Training the RF model for Tn 48h ##
registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        #repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_CN_Tn_48h <- train(CN~.,
                      data=CN_48h,
                      method="rf",
                      ntree=5000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(CN_48h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 21.57055 hours
print(rf_CN_Tn_48h)
print(rf_CN_Tn_48h$finalModel)

# Type of random forest: regression
# Number of trees: 5000
# No. of variables tried at each split: 931
# Mean of squared residuals: 0.4919163
# % Var explained: 50.98

saveRDS(rf_CN_Tn_48h, "models/protein/rf_CN_Tn_48h_Lib_V2_11092021.RDS")

gc()
rm(rf_CN_Tn_48h,CN_48h)


