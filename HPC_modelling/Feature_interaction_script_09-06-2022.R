library(plyr)
library(dplyr)
library(doMC)
library(biomaRt)
library(tidyverse)
library(caret)
library(e1071)
library(xgboost)




###__________________________________________________________________________________________###
###_______________________________________ immune subsets ___________________________________###
###__________________________________________________________________________________________###


setwd("~/Analysis/")

## Biomart ##
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version","ensembl_transcript_id","ensembl_gene_id","ccds","transcript_biotype"), mart = ensembl)
#
#
# ## parameters ##
print("importing lib")
Sequence_parameters_RNA <- read.delim("~/Analysis/libraries/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB (copy).csv", sep=";", dec=",")
dim(Sequence_parameters_RNA)





###__________________________________________________________________________________________###
###_________________________________________ B naive ________________________________________###
###__________________________________________________________________________________________###
# 
## prep data ##
# print("preping data for Tn 0h")
# TPM_Bnaive <- read.delim("~/Analysis/data/immune_cells_RNA/B_naive_TPM_log10.csv",sep = ";",dec = ",")
# TPM_Bnaive <- merge(TPM_Bnaive,Sequence_parameters_RNA,by.x="ID",by.y="tx.id")
# 
# 
# print("integrating feature counts per gene")
# registerDoMC(10)
# TPM_Bnaive_param <- ddply(TPM_Bnaive,"ensembl_gene_id", numcolwise(mean), .parallel = T, .progress = T)
# TPM_Bnaive_param$TPM <- NULL
# 
# print("integrating TPM per gene")
# TPM_Bnaive$TPM <- 10^(TPM_Bnaive$TPM)
# TPM_Bnaive <- data.frame("ID"=TPM_Bnaive$ensembl_gene_id ,"TPM"=TPM_Bnaive$TPM)
# TPM_Bnaive <- ddply(TPM_Bnaive,"ID", numcolwise(sum), .parallel = T, .progress = T)
# TPM_Bnaive <- merge(TPM_Bnaive,TPM_Bnaive_param,by.x="ID",by.y="ensembl_gene_id", all.x=T)
# 
# 
# TPM_Bnaive$TPM <- log10(TPM_Bnaive$TPM)
# 
# TPM_Bnaive[is.na(TPM_Bnaive)]= 0
# TPM_Bnaive <- subset(TPM_Bnaive, TPM_Bnaive$TPM>=-1)
# print(dim(TPM_Bnaive)) # 13729  2682
# 
# rownames(TPM_Bnaive) <- TPM_Bnaive$ID
# TPM_Bnaive$ID <- NULL
# TPM_Bnaive$ensembl_gene_id <- NULL
# TPM_Bnaive$gene_biotype <- NULL
# TPM_Bnaive$gene_name <- NULL
# 
# ## Feature selection ##
# 
# print("removing useless columns")
# TPM_Bnaive[is.na(TPM_Bnaive)]=0
# TPM_Bnaive <- TPM_Bnaive[, -caret::nearZeroVar(TPM_Bnaive, allowParallel = TRUE, uniqueCut = 0.1)]
# print(dim(TPM_Bnaive)) # 13729  2637
# 
# gc()
# 
# 
# ## test / train sets ##
# train_row <- sample(1:nrow(TPM_Bnaive), 0.8 * nrow(TPM_Bnaive))
# TPM_Bnaive_train <- TPM_Bnaive[train_row,]
# TPM_Bnaive_test <- TPM_Bnaive[-train_row,]
# 
# dim(TPM_Bnaive_train)
# dim(TPM_Bnaive_test)
# 
# write.table(TPM_Bnaive_train,"~/Analysis/split_dataset/TPM_Bnaive_train",sep=";")
# write.table(TPM_Bnaive_test,"~/Analysis/split_dataset/TPM_Bnaive_test",sep=";")

TPM_Bnaive_train <- read.delim("~/Analysis/split_dataset/TPM_Bnaive_train",sep=";")
TPM_Bnaive_test <- read.delim("~/Analysis/split_dataset/TPM_Bnaive_test",sep=";")

TPM_Bnaive_train_Else <- TPM_Bnaive_train %>% dplyr::select(-contains("UTR5"))
TPM_Bnaive_train_UTR5 <- TPM_Bnaive_train %>% dplyr::select(contains("UTR5"))

TPM_Bnaive_train_UTR5 <- TPM_Bnaive_train_UTR5[sample(1:nrow(TPM_Bnaive_train_UTR5)),]

TPM_Bnaive_train <- cbind.data.frame(TPM_Bnaive_train_UTR5,TPM_Bnaive_train_Else)





## Training the XGB model ##
print("modeling")
registerDoMC(1)
set.seed(12345)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        seeds = NULL) # NULL will set the seeds using a random set of integers (from doc)

xgbGrid <- expand.grid(nrounds = 10000,
                             max_depth = 6,
                             colsample_bytree = 0.4,
                             eta = 0.05,
                             gamma=0,
                             min_child_weight = 0.9,
                             subsample = 1)

start_time <- Sys.time()
xgb_TPM_Bnaive <- train(TPM~.,
                        data=TPM_Bnaive_train,
                        method="xgbTree",
                        trControl=control,
                        #metric="Rsquared",
                        tuneGrid= xgbGrid,
                        na.action = na.omit,
                        nthread=80,
                        verbose = TRUE)

end_time <- Sys.time()
print(end_time - start_time) # Time difference of  hours

print(xgb_TPM_Bnaive)

saveRDS(xgb_TPM_Bnaive, "~/Analysis/models/xgb_TPM_Bnaive_UTR5_dropped_03_05_2022.RDS")

# xgb_TPM_Tn_0h <- readRDS( "~/Analysis/models/xgb_TPM_Bnaive_03_05_2022.RDS")

xgb_TPM_Bnaive_predict_test <- predict(xgb_TPM_Bnaive,TPM_Bnaive_test)
xgb_TPM_Bnaive_predict_test <- data.frame(xgb_TPM_Bnaive_predict_test)

print(cor(xgb_TPM_Bnaive_predict_test$xgb_TPM_Bnaive_predict_test,TPM_Bnaive_test$TPM, method = "pearson", use = "complete.obs")^2)


gc()
rm(TPM_Bnaive,TPM_Bnaive_param)


xgb_TPM_Bnaive <- readRDS( "~/Analysis/models/xgb_TPM_Bnaive_03_05_2022.RDS")
xgb_TPM_Bnaive_no_5UTR <- readRDS( "~/Analysis/models/xgb_TPM_Bnaive_UTR5_dropped_03_05_2022.RDS")


var_imp_xgb_TPM_Bnaive <- data.frame(caret::varImp(xgb_TPM_Bnaive,scale=F)$importance)
var_imp_xgb_TPM_Bnaive_no_5UTR <- data.frame(caret::varImp(xgb_TPM_Bnaive_no_5UTR,scale=F)$importance)


var_imp_xgb_TPM_Bnaive$ID <- rownames(var_imp_xgb_TPM_Bnaive)
var_imp_xgb_TPM_Bnaive_no_5UTR$ID <- rownames(var_imp_xgb_TPM_Bnaive_no_5UTR)


var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR <- merge(var_imp_xgb_TPM_Bnaive,var_imp_xgb_TPM_Bnaive_no_5UTR,by="ID")
colnames(var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR) <- c("ID","var_imp_full","var_imp_no_5UTR")


var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_full <- log10(var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_full)
var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_no_5UTR <- log10(var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_no_5UTR)
var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_full <- as.numeric(gsub("-Inf","",var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_full))
var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_no_5UTR <- as.numeric(gsub("-Inf","",var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_no_5UTR))


var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$LFC_no_5UTR_vs_full <- var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_full-var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$var_imp_no_5UTR

var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$LFC_no_5UTR_vs_full


ggplot(var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR,aes(x=(var_imp_full),y=(var_imp_no_5UTR)))+
  geom_point()+
  scale_x_continuous(limits = c(-4,-1))+
  scale_y_continuous(limits = c(-4,-1))+
  geom_text(data=var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR[(var_imp_xgb_TPM_Bnaive_full_vs_no_5UTR$LFC_no_5UTR_vs_full)< -0.5,],aes(label=ID),nudge_y = .5)+
  theme(aspect.ratio = 1)



