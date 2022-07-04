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

## Biomart ##
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")

#All_filters <- listAttributes(ensembl)
#View(All_filters)
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)


## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###____________________________________________Tn_0h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_0h <- read.delim("~/Analysis/data/RNA/Tnaive_0h/quant.sf")

TPM_0h$Name <- mapply(strsplit(as.character(TPM_0h$Name),"\\."),FUN=function(x){(as.character(x)[1])})

TPM_0h <- data.frame("ID"=TPM_0h$Name,"TPM"=TPM_0h$TPM)
TPM_0h <- subset(TPM_0h, TPM_0h$TPM>0)
TPM_0h <- data.frame("ID"=TPM_0h$ID,"TPM"=TPM_0h$TPM)
TPM_0h <- merge(TPM_0h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")


print("integrating feature counts per gene")
registerDoMC(10)
TPM_0h_param <- ddply(TPM_0h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_0h_param$TPM <- NULL

print("integrating TPM per gene")
TPM_0h <- data.frame("ID"=TPM_0h$ensembl_gene_id ,"TPM"=TPM_0h$TPM)
TPM_0h <- ddply(TPM_0h,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_0h <- merge(TPM_0h,TPM_0h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_0h$TPM <- log10(TPM_0h$TPM)

TPM_0h[is.na(TPM_0h)]= 0
TPM_0h <- subset(TPM_0h, TPM_0h$TPM>=-1)
print(dim(TPM_0h)) # 13729  2682

rownames(TPM_0h) <- TPM_0h$ID
TPM_0h$ID <- NULL
TPM_0h$ensembl_gene_id <- NULL
TPM_0h$gene_biotype <- NULL
TPM_0h$gene_name <- NULL


write.table(TPM_0h,"~/Analysis/counts&libs/TPM_0h_with_lib_per_gene.csv",row.names = T)

## Feature selection ##
## Tn 0h ##

print("removing useless columns")
TPM_0h[is.na(TPM_0h)]=0
TPM_0h <- TPM_0h[, -caret::nearZeroVar(TPM_0h, allowParallel = TRUE, uniqueCut = dim(TPM_0h)[1]*0.01)]
print(dim(TPM_0h)) # 13729  2637

gc()


## Training the RF model for Tn 0h ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_Tn_0h <- train(TPM~.,
                      data=TPM_0h,
                      method="rf",
                      ntree=3000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_0h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours

print(rf_TPM_Tn_0h)

print(rf_TPM_Tn_0h$finalModel)

saveRDS(rf_TPM_Tn_0h, "~/Analysis/models/rf_deep_3000trees_TPM_Tn_0h_11_09_2021.RDS")

gc()
rm(TPM_0h,TPM_0h_param)





###__________________________________________________________________________________________###
###____________________________________________Tn_6h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_6h <- read.delim("~/Analysis/data/RNA/Tnaive_act6h/quant.sf")
TPM_6h$Name <- mapply(strsplit(as.character(TPM_6h$Name),"\\."),FUN=function(x){(as.character(x)[1])})

TPM_6h <- data.frame("ID"=TPM_6h$Name,"TPM"=TPM_6h$TPM)
TPM_6h <- subset(TPM_6h, TPM_6h$TPM>0)
TPM_6h <- data.frame("ID"=TPM_6h$ID,"TPM"=TPM_6h$TPM)
TPM_6h <- merge(TPM_6h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")

print("integrating feature counts per gene")
registerDoMC(10)
TPM_6h_param <- ddply(TPM_6h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_6h_param$TPM <- NULL

print("integrating TPM per gene")
TPM_6h <- data.frame("ID"=TPM_6h$ensembl_gene_id ,"TPM"=TPM_6h$TPM)
TPM_6h <- ddply(TPM_6h,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_6h <- merge(TPM_6h,TPM_6h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)

TPM_6h$TPM <- log10(TPM_6h$TPM)
TPM_6h[is.na(TPM_6h)]= 0
TPM_6h <- subset(TPM_6h, TPM_6h$TPM>=-1)
print(dim(TPM_6h)) # 13729  2682

rownames(TPM_6h) <- TPM_6h$ID
TPM_6h$ID <- NULL
TPM_6h$ensembl_gene_id <- NULL
TPM_6h$gene_biotype <- NULL
TPM_6h$gene_name <- NULL


write.table(TPM_6h,"~/Analysis/counts&libs/TPM_6h_with_lib_per_gene.csv",row.names = T)


## Feature selection ##
## Tn 6h ##
print("removing useless columns")
TPM_6h[is.na(TPM_6h)]=0
TPM_6h <- TPM_6h[, -caret::nearZeroVar(TPM_6h, allowParallel = TRUE, uniqueCut = dim(TPM_6h)[1]*0.01)]
print(dim(TPM_6h)) # 13729  2637

gc()

## Training the RF model for Tn 6h ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_Tn_6h <- train(TPM~.,
                      data=TPM_6h,
                      method="rf",
                      ntree=3000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_6h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours
print(rf_TPM_Tn_6h)
print(rf_TPM_Tn_6h$finalModel)
saveRDS(rf_TPM_Tn_6h, "~/Analysis/models/rf_deep_3000trees_TPM_Tn_6h_11_09_2021.RDS")


gc()
rm(TPM_6h,TPM_6h_param)




###__________________________________________________________________________________________###
###____________________________________________Tn_24h_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_24h <- read.delim("~/Analysis/data/RNA/Tnaive_act24h/quant.sf")
TPM_24h$Name <- mapply(strsplit(as.character(TPM_24h$Name),"\\."),FUN=function(x){(as.character(x)[1])})

TPM_24h <- data.frame("ID"=TPM_24h$Name,"TPM"=TPM_24h$TPM)
TPM_24h <- subset(TPM_24h, TPM_24h$TPM>0)
TPM_24h <- data.frame("ID"=TPM_24h$ID,"TPM"=TPM_24h$TPM)
TPM_24h <- merge(TPM_24h,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")

print("integrating feature counts per gene")
registerDoMC(10)
TPM_24h_param <- ddply(TPM_24h,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_24h_param$TPM <- NULL

print("integrating TPM per gene")
TPM_24h <- data.frame("ID"=TPM_24h$ensembl_gene_id ,"TPM"=TPM_24h$TPM)
TPM_24h <- ddply(TPM_24h,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_24h <- merge(TPM_24h,TPM_24h_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
TPM_24h$TPM <- log10(TPM_24h$TPM)

TPM_24h[is.na(TPM_24h)]= 0
TPM_24h <- subset(TPM_24h, TPM_24h$TPM>=-1)
print(dim(TPM_24h)) # 13729  2682

rownames(TPM_24h) <- TPM_24h$ID
TPM_24h$ID <- NULL
TPM_24h$ensembl_gene_id <- NULL
TPM_24h$gene_biotype <- NULL
TPM_24h$gene_name <- NULL


write.table(TPM_24h,"~/Analysis/counts&libs/TPM_24h_with_lib_per_gene.csv",row.names = T)


## Feature selection ##
## Tn 24h ##

print("removing useless columns")
TPM_24h[is.na(TPM_24h)]=0
TPM_24h <- TPM_24h[, -caret::nearZeroVar(TPM_24h, allowParallel = TRUE, uniqueCut = dim(TPM_24h)[1]*0.01)]
print(dim(TPM_24h)) # 13729  2637

gc()

## Training the RF model for Tn 24h ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_Tn_24h <- train(TPM~.,
                      data=TPM_24h,
                      method="rf",
                      ntree=3000,
                      trControl=control,
                      metric="Rsquared",
                      tuneGrid= data.frame(mtry = round((ncol(TPM_24h)-1)/3)),
                      na.action = na.omit,
                      importance = TRUE,
                      verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours
print(rf_TPM_Tn_24h)
print(rf_TPM_Tn_24h$finalModel)
saveRDS(rf_TPM_Tn_24h, "~/Analysis/models/rf_deep_3000trees_TPM_Tn_24h_11_09_2021.RDS")

gc()
rm(TPM_24h,TPM_24h_param)



rf_TPM_Tn_24h <- readRDS("~/Analysis/models/rf_deep_3000trees_TPM_Tn_24h_11_09_2021.RDS")


rf_TPM_Tn_24h_IMP <- varImp(rf_TPM_Tn_24h$finalModel)
rf_TPM_Tn_24h_IMP <- data.frame(rf_TPM_Tn_24h_IMP)
rf_TPM_Tn_24h_IMP$ID <- rownames(rf_TPM_Tn_24h_IMP) 

ggplot(TPM_24h, aes(x=total_m7G,y=TPM))+
  geom_point(alpha=0.1,stroke=0)+
  geom_smooth(method="loess")



###__________________________________________________________________________________________###
###____________________________________________CD8_Tn_________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_CD8_Tn <- read.delim("~/Analysis/data/RNA/CD8_Tn_TPM_log10.csv",sep=";",dec=",")
TPM_CD8_Tn <- merge(TPM_CD8_Tn,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")

print("integrating feature counts per gene")
registerDoMC(10)
TPM_CD8_Tn_param <- ddply(TPM_CD8_Tn,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_CD8_Tn_param$TPM <- NULL

print("integrating TPM per gene")
TPM_CD8_Tn <- data.frame("ID"=TPM_CD8_Tn$ensembl_gene_id ,"TPM"=TPM_CD8_Tn$TPM)
TPM_CD8_Tn <- ddply(TPM_CD8_Tn,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Tn <- merge(TPM_CD8_Tn,TPM_CD8_Tn_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_CD8_Tn$TPM <- log10(TPM_CD8_Tn$TPM)

TPM_CD8_Tn[is.na(TPM_CD8_Tn)]= 0
TPM_CD8_Tn <- subset(TPM_CD8_Tn, TPM_CD8_Tn$TPM>=-1)
print(dim(TPM_CD8_Tn)) # 13729  2682

rownames(TPM_CD8_Tn) <- TPM_CD8_Tn$ID
TPM_CD8_Tn$ID <- NULL
TPM_CD8_Tn$ensembl_gene_id <- NULL
TPM_CD8_Tn$gene_biotype <- NULL
TPM_CD8_Tn$gene_name <- NULL

write.table(TPM_CD8_Tn,"~/Analysis/counts&libs/TPM_CD8_Tn_with_lib_per_gene.csv",row.names = T)


## Feature selection ##
## Tn CD8_Tn ##

print("removing useless columns")
TPM_CD8_Tn[is.na(TPM_CD8_Tn)]=0
TPM_CD8_Tn_TPM <- TPM_CD8_Tn$TPM
TPM_CD8_Tn <- TPM_CD8_Tn[, -caret::nearZeroVar(TPM_CD8_Tn, allowParallel = TRUE, uniqueCut = dim(TPM_CD8_Tn)[1]*0.01)]
TPM_CD8_Tn$TPM <- TPM_CD8_Tn_TPM
print(dim(TPM_CD8_Tn)) # 13729  2637

gc()


## Training the RF model for Tn CD8_Tn ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_CD8_Tn <- train(TPM~.,
                       data=TPM_CD8_Tn,
                       method="rf",
                       ntree=3000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(TPM_CD8_Tn)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours
print(rf_TPM_CD8_Tn)
print(rf_TPM_CD8_Tn$finalModel)
saveRDS(rf_TPM_CD8_Tn, "~/Analysis/models/rf_deep_3000trees_TPM_CD8_Tn_11_09_2021.RDS")

gc()
rm(TPM_CD8_Tn,TPM_CD8_Tn_param)




###__________________________________________________________________________________________###
###____________________________________________CD8_Tcm________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_CD8_Tcm <- read.delim("~/Analysis/data/RNA/CD8_Tcm_TPM_log10.csv",sep=";",dec=",")
TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
dim(TPM_CD8_Tcm)

print("integrating feature counts per gene")
registerDoMC(10)
TPM_CD8_Tcm_param <- ddply(TPM_CD8_Tcm,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_CD8_Tcm_param$TPM <- NULL

print("integrating TPM per gene")
TPM_CD8_Tcm <- data.frame("ID"=TPM_CD8_Tcm$ensembl_gene_id ,"TPM"=TPM_CD8_Tcm$TPM)
TPM_CD8_Tcm <- ddply(TPM_CD8_Tcm,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Tcm <- merge(TPM_CD8_Tcm,TPM_CD8_Tcm_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_CD8_Tcm$TPM <- log10(TPM_CD8_Tcm$TPM)
TPM_CD8_Tcm[is.na(TPM_CD8_Tcm)]= 0
TPM_CD8_Tcm <- subset(TPM_CD8_Tcm, TPM_CD8_Tcm$TPM>=-1)
print(dim(TPM_CD8_Tcm)) # 13729  2682

rownames(TPM_CD8_Tcm) <- TPM_CD8_Tcm$ID
TPM_CD8_Tcm$ID <- NULL
TPM_CD8_Tcm$ensembl_gene_id <- NULL
TPM_CD8_Tcm$gene_biotype <- NULL
TPM_CD8_Tcm$gene_name <- NULL

write.table(TPM_CD8_Tcm,"~/Analysis/counts&libs/TPM_CD8_Tcm_with_lib_per_gene.csv",row.names = T)

## Feature selection ##
## Tn CD8_Tcm ##

print("removing useless columns")
TPM_CD8_Tcm[is.na(TPM_CD8_Tcm)]=0
TPM_CD8_Tcm_TPM <- TPM_CD8_Tcm$TPM
TPM_CD8_Tcm <- TPM_CD8_Tcm[, -caret::nearZeroVar(TPM_CD8_Tcm, allowParallel = TRUE, uniqueCut = dim(TPM_CD8_Tcm)[1]*0.01)]
TPM_CD8_Tcm$TPM <- TPM_CD8_Tcm$TPM
print(dim(TPM_CD8_Tcm)) # 13729  2637

gc()

## Training the RF model for Tn CD8_Tcm ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_CD8_Tcm <- train(TPM~.,
                       data=TPM_CD8_Tcm,
                       method="rf",
                       ntree=3000,
                       trControl=control,
                       metric="Rsquared",
                       tuneGrid= data.frame(mtry = round((ncol(TPM_CD8_Tcm)-1)/3)),
                       na.action = na.omit,
                       importance = TRUE,
                       verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours
print(rf_TPM_CD8_Tcm)
print(rf_TPM_CD8_Tcm$finalModel)
saveRDS(rf_TPM_CD8_Tcm, "~/Analysis/models/rf_deep_3000trees_TPM_CD8_Tcm_11_09_2021.RDS")


gc()
rm(TPM_CD8_Tcm,TPM_CD8_Tcm_param)



###__________________________________________________________________________________________###
###____________________________________________CD8_Tem________________________________________###
###__________________________________________________________________________________________###

## prep data ##
print("preping data")
TPM_CD8_Tem <- read.delim("~/Analysis/data/RNA/CD8_Tem_TPM_log10.csv",sep=";",dec=",")

#dim(TPM_CD8_Tem)
TPM_CD8_Tem <- merge(TPM_CD8_Tem,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
dim(TPM_CD8_Tem)

print("integrating feature counts per gene")
registerDoMC(10)
TPM_CD8_Tem_param <- ddply(TPM_CD8_Tem,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
TPM_CD8_Tem_param$TPM <- NULL

print("integrating TPM per gene")
TPM_CD8_Tem <- data.frame("ID"=TPM_CD8_Tem$ensembl_gene_id ,"TPM"=TPM_CD8_Tem$TPM)
TPM_CD8_Tem <- ddply(TPM_CD8_Tem,"ID", numcolwise(sum), .parallel = T, .progress = T)
TPM_CD8_Tem <- merge(TPM_CD8_Tem,TPM_CD8_Tem_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)


TPM_CD8_Tem$TPM <- log10(TPM_CD8_Tem$TPM)

TPM_CD8_Tem[is.na(TPM_CD8_Tem)]= 0
TPM_CD8_Tem <- subset(TPM_CD8_Tem, TPM_CD8_Tem$TPM>=-1)
print(dim(TPM_CD8_Tem)) # 13729  2682

rownames(TPM_CD8_Tem) <- TPM_CD8_Tem$ID
TPM_CD8_Tem$ID <- NULL
TPM_CD8_Tem$ensembl_gene_id <- NULL
TPM_CD8_Tem$gene_biotype <- NULL
TPM_CD8_Tem$gene_name <- NULL

write.table(TPM_CD8_Tem,"~/Analysis/counts&libs/TPM_CD8_Tem_with_lib_per_gene.csv",row.names = T)

## Feature selection ##
## Tn CD8_Tem ##

print("removing useless columns")
TPM_CD8_Tem[is.na(TPM_CD8_Tem)]=0
TPM_CD8_Tem_TPM <- TPM_CD8_Tem$TPM
TPM_CD8_Tem <- TPM_CD8_Tem[, -caret::nearZeroVar(TPM_CD8_Tem, allowParallel = TRUE, uniqueCut = dim(TPM_CD8_Tem)[1]*0.01)]
TPM_CD8_Tem$TPM <- TPM_CD8_Tem_TPM
print(dim(TPM_CD8_Tem)) # 13729  2637

gc()

## Training the RF model for Tn CD8_Tem ##
print("modeling")

registerDoMC(10)
set.seed(12345)
control <- trainControl(method="cv",
                        number=10,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)


start_time <- Sys.time()
rf_TPM_CD8_Tem <- train(TPM~.,
                        data=TPM_CD8_Tem,
                        method="rf",
                        ntree=3000,
                        trControl=control,
                        metric="Rsquared",
                        tuneGrid= data.frame(mtry = round((ncol(TPM_CD8_Tem)-1)/3)),
                        na.action = na.omit,
                        importance = TRUE,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.78211 hours
print(rf_TPM_CD8_Tem)
print(rf_TPM_CD8_Tem$finalModel)
saveRDS(rf_TPM_CD8_Tem, "~/Analysis/models/rf_deep_3000trees_TPM_CD8_Tem_11_09_2021.RDS")


gc()
rm(TPM_CD8_Tem,TPM_CD8_Tem_param)


