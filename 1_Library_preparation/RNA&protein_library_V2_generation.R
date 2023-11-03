#title: "Sequence feature library generation"
#author: "Benoit P Nicolet"
#date: "02/09/2023"


###__________________________________________________________________________________________###
###______________________________________Initialization______________________________________###
###__________________________________________________________________________________________###

library(plyr)
library(dplyr)
library(seqinr)
library(stringr)
library(qdapRegex)
library(dplyr)
library(Biostrings)
library(biomaRt)
# library(ggplot2)
library(stringdist)
library(fuzzyjoin)
library(coRdon)


setwd("~/") # Setting up the work directory 


## Acquiring Ensembl annotations ##

ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
IDs_genenames_coding <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","transcript_biotype"), mart = ensembl)
IDs_genenames_coding <- IDs_genenames_coding[IDs_genenames_coding$transcript_biotype=="protein_coding",]


GeneID_uniprot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","uniprotswissprot","transcript_biotype"), mart = ensembl)
GeneID_uniprot <- GeneID_uniprot[GeneID_uniprot$transcript_biotype=="protein_coding",]
GeneID_uniprot$transcript_biotype <- NULL
# dim(GeneID_uniprot)




## Importing ATTRACT DB

## ATTRACT db downloaded on June 11th 2020 from https://attract.cnic.es/index
ATTRACT <- read.delim("./ATtRACT_db.txt")
ATTRACT <- subset(ATTRACT,ATTRACT$Organism=="Homo_sapiens")
dim(ATTRACT) # 3256 motifs

## converting Us in Ts
ATTRACT$Motif <- gsub("U","T",ATTRACT$Motif)

table(duplicated(ATTRACT$Motif))
ATTRACT <- ATTRACT[!duplicated(ATTRACT$Motif),] ## Removing duplicates motifs
dim(ATTRACT) # 2271 motifs

dim(ATTRACT[!duplicated(ATTRACT$Gene_name),]) # Should be == dim(ATTRACT)




###__________________________________________________________________________________________###
###__________________________________________3'UTR___________________________________________###
###__________________________________________________________________________________________###


# Here we re-format the fasta format to table format. Run once ##
UTR3_fasta <- readtext::readtext("./3UTR_ensembl_r104.txt")
glimpse(UTR3_fasta)
UTR3_fasta <- gsub("\\\n","",UTR3_fasta) # removing the newlines
UTR3_fasta <- gsub(">","\\\n>",UTR3_fasta) ## replace > by \n>

write(UTR3_fasta,"./3UTR_ensembl_r104_formatted.txt")


## importing the formatted Fasta file and processing ##

UTR3_fasta <- read.delim("./3UTR_ensembl_r104_formatted.txt",sep="\t",header = F)
UTR3_fasta$sequence <- UTR3_fasta$V1
colnames(UTR3_fasta)[1] <- "tx.id"
UTR3_fasta$sequence <- gsub('^................','',UTR3_fasta$sequence)


UTR3_fasta$tx.id <- gsub(">","",UTR3_fasta$tx.id)
UTR3_fasta$tx.id <- gsub("A","",UTR3_fasta$tx.id)
UTR3_fasta$tx.id <- gsub("T","",UTR3_fasta$tx.id)
UTR3_fasta$tx.id <- gsub("G","",UTR3_fasta$tx.id)
UTR3_fasta$tx.id <- gsub("C","",UTR3_fasta$tx.id)
UTR3_fasta$tx.id <- gsub("ENS","ENST",UTR3_fasta$tx.id)

UTR3_fasta$UTR3_length <- mapply(strsplit(as.character(UTR3_fasta$sequence),","),FUN=function(x){nchar(x)}) # creating a column with 3UTR length
# dim(UTR3_fasta)

UTR3_fasta <- UTR3_fasta[!UTR3_fasta$sequence=="Sequence unavailable",]
table(duplicated(UTR3_fasta$tx.id)) # should be 0







## Counting motifs ##

# Here I made a for loop to count the motifs contained in ATTRACT dataframe in the UTR3_fasta sequence #
# I also implemented a rough progression counter per 10% #

datalist = list() # creating the list that will contain the data

{
  start_time <- Sys.time()
  for (i in 1:length(ATTRACT$Motif)) {
    dat <- as.vector(str_count(UTR3_fasta$sequence, pattern=as.character(ATTRACT$Motif[i])))
    datalist[[i]] <- dat # add it to your list
    
    ## This next part is just to return the progress (per 10% increase)
    k <- round((as.numeric(i-1)/length(ATTRACT$Motif))*100,-1)
    if (round((as.numeric(i)/length(ATTRACT$Motif))*100,-1)>k)
    {print((paste0(round((as.numeric(i)/length(ATTRACT$Motif))*100,-1),"%")))}
    else {}
    flush.console()
  }

  end_time <- Sys.time()
  end_time - start_time # takes ~10-30min 
}

UTR3_database <- do.call(cbind,datalist) # transforming the list into a df

colnames(UTR3_database) <- paste0(ATTRACT$Motif) # defining motifs as column names
rownames(UTR3_database) <- UTR3_fasta$tx.id # plugging the tx id back
# dim(UTR3_database) # QC


UTR3_database <- UTR3_database[,order(colnames(UTR3_database))] # ordering the columns for convienience 
UTR3_database_export <- data.frame(UTR3_database) # creating a new object for safety
UTR3_database_export$tx.id <- rownames(UTR3_database_export)


# Not all motifs are found in the transcriptome, so I remove the empty columns #
dim(UTR3_database[, colSums(UTR3_database) !=0 ]) # checking: 2235 non-empty columns 
UTR3_database <- UTR3_database[, colSums(UTR3_database) != 0] # removing empty columns
UTR3_database <- UTR3_database[,order(colnames(UTR3_database))] # reordering

UTR3_fasta <- cbind(UTR3_fasta,UTR3_database) # combining the 3UTR sequence, length, and sequence motifs occurence


## including G-C content ##

UTR3_fasta$G_count <- str_count(UTR3_fasta$sequence,pattern ="G") # counting Gs
UTR3_fasta$C_count <- str_count(UTR3_fasta$sequence,pattern ="C") # counting Cs
UTR3_fasta$GC_percentage <- ((UTR3_fasta$G_count + UTR3_fasta$C_count)/UTR3_fasta$UTR3_length)*100 # making a G+C %

# droping columns
UTR3_fasta$G_count <- NULL
UTR3_fasta$C_count <- NULL



## making tables pretty for export ##
UTR3_fasta_export <- UTR3_fasta
UTR3_fasta_export$sequence <- NULL # removing the sequence for export
colnames(UTR3_fasta_export) <- paste0(colnames(UTR3_fasta_export),"_UTR3") # To make the column names unique, I just append "_UTR3" to the original column names.

# Now we export the 3'UTR info to file for future use #
write.table(UTR3_fasta_export,"./3UTR_library_v2_RBP_GC_length.csv",sep = ";", dec = ",",row.names = F)




###__________________________________________________________________________________________###
###____________________________________________CDS___________________________________________###
###__________________________________________________________________________________________###

# Here I run the same script than 3'UTR for CDS region, refer to the annotation of 3'UTR #

## importing the formated Fasta file and processing ##

## Here we re-format the fasta format to table format ##
# Run once # 
CDS_fasta <- readtext::readtext("./CDS_ensembl_r104.txt")
glimpse(CDS_fasta)
CDS_fasta <- gsub("\\\n","",CDS_fasta) # removing the newlines
CDS_fasta <- gsub(">","\\\n>",CDS_fasta) ## replace > by \n>
write(CDS_fasta,"./CDS_ensembl_r104_formatted.txt")



CDS_fasta <- read.delim("./CDS_ensembl_r104_formatted.txt",sep="\t",header = F)


CDS_fasta$sequence <- CDS_fasta$V1
colnames(CDS_fasta)[1] <- "tx.id"
CDS_fasta$sequence <- gsub('^................','',CDS_fasta$sequence)


CDS_fasta$tx.id <- gsub(">","",CDS_fasta$tx.id)
CDS_fasta$tx.id <- gsub("A","",CDS_fasta$tx.id)
CDS_fasta$tx.id <- gsub("T","",CDS_fasta$tx.id)
CDS_fasta$tx.id <- gsub("G","",CDS_fasta$tx.id)
CDS_fasta$tx.id <- gsub("C","",CDS_fasta$tx.id)
CDS_fasta$tx.id <- gsub("ENS","ENST",CDS_fasta$tx.id)

CDS_fasta$CDS_length <- mapply(strsplit(as.character(CDS_fasta$sequence),","),FUN=function(x){nchar(x)})
CDS_fasta <- CDS_fasta[!CDS_fasta$sequence=="Sequence unavailable",]
table(duplicated(CDS_fasta$tx.id))




## Counting motifs ##
datalist_CDS = list()

{
  start_time <- Sys.time()
  for (i in 1:length(ATTRACT$Motif)) {
    dat <- as.vector(str_count(CDS_fasta$sequence, pattern=as.character(ATTRACT$Motif[i])))
    datalist_CDS[[i]] <- dat # add it to your list
    ## This next part is just to return the progress (per 10% increase)
    k <- round((as.numeric(i-1)/length(ATTRACT$Motif))*100,-1)
    if (round((as.numeric(i)/length(ATTRACT$Motif))*100,-1)>k)
    {print((paste0(round((as.numeric(i)/length(ATTRACT$Motif))*100,-1),"%")))}
    else {}
    flush.console()
  }
  
  end_time <- Sys.time()
  end_time - start_time
}

CDS_database <- do.call(cbind,datalist_CDS)
colnames(CDS_database) <- paste0(ATTRACT$Motif)
rownames(CDS_database) <- CDS_fasta$tx.id
dim(CDS_database)


CDS_database <- CDS_database[,order(colnames(CDS_database))]
CDS_database_export <- data.frame(CDS_database)
CDS_database_export$tx.id <- CDS_fasta$tx.id
dim(CDS_database[, colSums(CDS_database) !=0 ]) # 2268 non-empty columns
CDS_database <- CDS_database[, colSums(CDS_database) != 0]
CDS_database <- CDS_database[,order(colnames(CDS_database))]
CDS_fasta <- cbind(CDS_fasta,CDS_database)
print(dim(CDS_fasta))

## G-C content ##
CDS_fasta$G_count <- str_count(CDS_fasta$sequence,pattern ="G")
CDS_fasta$C_count <- str_count(CDS_fasta$sequence,pattern ="C")
CDS_fasta$GC_percentage <- ((CDS_fasta$G_count + CDS_fasta$C_count)/CDS_fasta$CDS_length)*100
CDS_fasta$G_count <- NULL
CDS_fasta$C_count <- NULL




## Codon usage ##

## Getting all codon usage info! ##
CDS_fasta_codon <- data.frame(row.names = CDS_fasta$tx.id,"seq"=CDS_fasta$sequence) # making a data frame for DNAStringSet 
CDS_fasta_codon_count <- Biostrings::DNAStringSet(CDS_fasta_codon$seq,use.names = T) # making a DNAStringSet object
CDS_fasta_codon_count_table <- codonTable(CDS_fasta_codon_count) # here I get the codon usage of the DNAStringSet object
CDS_fasta_codon_count <- data.frame(codonCounts(CDS_fasta_codon_count_table)) # making a df of the data

CDS_fasta_codon_count$AA_length <- rowSums(CDS_fasta_codon_count) #computing the amino acid length
CDS_fasta_codon_count[1:64] <- (CDS_fasta_codon_count[1:64]/CDS_fasta_codon_count$AA_length)*100 # computing the frequency per codon
CDS_fasta_codon_count$AA_length <- NULL # dropping the AA length

colnames(CDS_fasta_codon_count) <- paste0(colnames(CDS_fasta_codon_count),"_codon") # to differentiate codons from sequence motifs, I add "_codon" to column names
CDS_fasta <- cbind.data.frame(CDS_fasta,CDS_fasta_codon_count) # adding the data to the CDS info



### Amino Acid usage ###
CDS_fasta$sequence <- gsub("N","A",CDS_fasta$sequence) #removing N characters
CDS_fasta$translation <- mapply(strsplit(as.character(CDS_fasta$sequence),","),FUN=function(x){as.character(translate(Biostrings::DNAString(x)))}) #translation of the CDS 

# here I use a crude but straight forward way to count amino acids #
CDS_fasta$STP <- str_count(CDS_fasta$translation,pattern ="\\*")
CDS_fasta$A_amino <- str_count(CDS_fasta$translation,pattern ="A")
CDS_fasta$C_amino <- str_count(CDS_fasta$translation,pattern ="C")
CDS_fasta$D_amino <- str_count(CDS_fasta$translation,pattern ="D")
CDS_fasta$E_amino <- str_count(CDS_fasta$translation,pattern ="E")
CDS_fasta$F_amino <- str_count(CDS_fasta$translation,pattern ="F")
CDS_fasta$G_amino <- str_count(CDS_fasta$translation,pattern ="G")
CDS_fasta$H_amino <- str_count(CDS_fasta$translation,pattern ="H")
CDS_fasta$I_amino <- str_count(CDS_fasta$translation,pattern ="I")
CDS_fasta$K_amino <- str_count(CDS_fasta$translation,pattern ="K")
CDS_fasta$L_amino <- str_count(CDS_fasta$translation,pattern ="L")
CDS_fasta$M_amino <- str_count(CDS_fasta$translation,pattern ="M")
CDS_fasta$N_amino <- str_count(CDS_fasta$translation,pattern ="N")
CDS_fasta$P_amino <- str_count(CDS_fasta$translation,pattern ="P")
CDS_fasta$Q_amino <- str_count(CDS_fasta$translation,pattern ="Q")
CDS_fasta$R_amino <- str_count(CDS_fasta$translation,pattern ="R")
CDS_fasta$S_amino <- str_count(CDS_fasta$translation,pattern ="S")
CDS_fasta$T_amino <- str_count(CDS_fasta$translation,pattern ="T")
CDS_fasta$V_amino <- str_count(CDS_fasta$translation,pattern ="V")
CDS_fasta$W_amino <- str_count(CDS_fasta$translation,pattern ="W")
CDS_fasta$Y_amino <- str_count(CDS_fasta$translation,pattern ="Y")

CDS_fasta$AA_length <- mapply(strsplit(as.character(CDS_fasta$translation),","),FUN=function(x){nchar(x)}) # CDS length in AA
CDS_fasta$STP <- NULL # I remove the stop codon info
CDS_fasta[2339:2358] <- (CDS_fasta[2339:2358]/CDS_fasta$AA_length)*100 # Computing the frequency of AA per length. The [2339:2358] is used to isolate the AA column and not the entire df



## making tables pretty for export}
CDS_fasta_export <- CDS_fasta
CDS_fasta_export$sequence <- NULL
colnames(CDS_fasta_export) <- paste0(colnames(CDS_fasta_export),"_CDS")

write.table(CDS_fasta_export,"./CDS_library_v2_RBP_GC_length.csv",sep = ";", dec = ",",row.names = F)




###__________________________________________________________________________________________###
###____________________________________________UTR5___________________________________________###
###__________________________________________________________________________________________###

## importing the formated Fasta file and processing ##

## Here we re-format the fasta format to table format ##
UTR5_fasta <- readtext::readtext("./5UTR_ensembl_r104.txt")
UTR5_fasta <- gsub("\\\n","",UTR5_fasta) # removing the newlines
UTR5_fasta <- gsub(">","\\\n>",UTR5_fasta) ## replace > by \n>
write(UTR5_fasta,"./5UTR_ensembl_r104_formatted.txt")

UTR5_fasta <- read.delim("./5UTR_ensembl_r104_formatted.txt",sep="\t",header = F)


UTR5_fasta$sequence <- UTR5_fasta$V1
colnames(UTR5_fasta)[1] <- "tx.id"
UTR5_fasta$sequence <- gsub('^................','',UTR5_fasta$sequence)
UTR5_fasta$tx.id <- gsub(">","",UTR5_fasta$tx.id)
UTR5_fasta$tx.id <- gsub("A","",UTR5_fasta$tx.id)
UTR5_fasta$tx.id <- gsub("T","",UTR5_fasta$tx.id)
UTR5_fasta$tx.id <- gsub("G","",UTR5_fasta$tx.id)
UTR5_fasta$tx.id <- gsub("C","",UTR5_fasta$tx.id)
UTR5_fasta$tx.id <- gsub("ENS","ENST",UTR5_fasta$tx.id)

UTR5_fasta$UTR5_length <- mapply(strsplit(as.character(UTR5_fasta$sequence),","),FUN=function(x){nchar(x)})
UTR5_fasta <- UTR5_fasta[!UTR5_fasta$sequence=="Sequence unavailable",]
table(duplicated(UTR5_fasta$tx.id))




## Counting motifs ##
datalist_UTR5 = list()

{
  start_time <- Sys.time()
  for (i in 1:length(ATTRACT$Motif)) {
    dat <- as.vector(str_count(UTR5_fasta$sequence, pattern=as.character(ATTRACT$Motif[i])))
    datalist_UTR5[[i]] <- dat # add it to your list
    ## This next part is just to return the progress (per 10% increase)
    k <- round((as.numeric(i-1)/length(ATTRACT$Motif))*100,-1)
    if (round((as.numeric(i)/length(ATTRACT$Motif))*100,-1)>k)
    {print((paste0(round((as.numeric(i)/length(ATTRACT$Motif))*100,-1),"%")))}
    else {}
    flush.console()
  }

  end_time <- Sys.time()
  end_time - start_time
}

UTR5_database <- do.call(cbind,datalist_UTR5)
colnames(UTR5_database) <- paste0(ATTRACT$Motif)
rownames(UTR5_database) <- UTR5_fasta$tx.id
dim(UTR5_database)

UTR5_database <- UTR5_database[,order(colnames(UTR5_database))]
UTR5_database_export <- data.frame(UTR5_database)
UTR5_database_export$tx.id <- UTR5_fasta$tx.id
dim(UTR5_database[, colSums(UTR5_database) !=0 ]) # 2226 non-empty columns
UTR5_database <- UTR5_database[, colSums(UTR5_database) != 0]
UTR5_database <- UTR5_database[,order(colnames(UTR5_database))]
UTR5_fasta <- cbind(UTR5_fasta,UTR5_database)


## G-C content ##
UTR5_fasta$G_count <- str_count(UTR5_fasta$sequence,pattern ="G")
UTR5_fasta$C_count <- str_count(UTR5_fasta$sequence,pattern ="C")
UTR5_fasta$GC_percentage <- ((UTR5_fasta$G_count + UTR5_fasta$C_count)/UTR5_fasta$UTR5_length)*100
UTR5_fasta$G_count <- NULL
UTR5_fasta$C_count <- NULL


## making tables pretty for export}
UTR5_fasta_export <- UTR5_fasta
UTR5_fasta_export$sequence <- NULL
colnames(UTR5_fasta_export) <- paste0(colnames(UTR5_fasta_export),"_UTR5")

write.table(UTR5_fasta_export,"./5UTR_library_v2_RBP_GC_length.csv",sep = ";", dec = ",",row.names = F)



###__________________________________________________________________________________________###
###____________________________________RNA modifications_____________________________________###
###__________________________________________________________________________________________###

## m6A content
# # Acquired from m6avar.renlab.org on SEP 02, 2021 # #
m6A <- read.delim("./RMVar_Human_basic_info_m6A.txt",sep="\t")
m6A <- data.frame("gene_name"=m6A$gene,"gene_id"=m6A$gene_database_id,"m6A_gene_region"=m6A$gene_region,"total_m6A"=1)
m6A <- m6A[!m6A$gene_name=="-",]
m6A <- m6A[grepl("Ensembl:",m6A$gene_id)==TRUE,]
m6A <- m6A[!grepl("intron",m6A$m6A_gene_region)==TRUE,]

# Here I bin the modifications per mRNA regions: #
m6A$CDS_m6a <- ifelse(grepl("CDS",m6A$m6A_gene_region)==TRUE,1,0)
m6A$UTR5_m6a <- ifelse(grepl("5'UTR",m6A$m6A_gene_region)==TRUE,1,0)
m6A$UTR3_m6a <- ifelse(grepl("3'UTR",m6A$m6A_gene_region)==TRUE,1,0)
m6A$gene_id <- mapply(strsplit(as.character(m6A$gene_id),","),FUN=function(x){(as.character(x)[1])})
m6A$gene_id <- gsub("Ensembl\\:","",m6A$gene_id)
m6A$m6A_gene_region <- NULL

m6A <- ddply(m6A,"gene_id", numcolwise(sum))
dim(m6A) # 38730 genes with m6A sites






### m5C ###
m5C <- read.delim("./RMVar_Human_basic_info_m5C.txt",sep="\t")
m5C <- m5C[!m5C$gene=="-",]
m5C <- m5C[grepl("Ensembl:",m5C$gene_database_id)==TRUE,]
m5C <- m5C[!grepl("intron",m5C$gene_region)==TRUE,]
m5C <- m5C[grepl("Protein coding",m5C$gene_type)==TRUE,]
m5C <- data.frame("gene_name"=m5C$gene,"gene_id"=m5C$gene_database_id,"gene_region"=m5C$gene_region,"total_m5C"=1)
dim(m5C)

m5C$CDS_m5C <- ifelse(grepl("CDS",m5C$gene_region)==TRUE,1,0)
m5C$UTR5_m5C <- ifelse(grepl("5'UTR",m5C$gene_region)==TRUE,1,0)
m5C$UTR3_m5C <- ifelse(grepl("3'UTR",m5C$gene_region)==TRUE,1,0)

m5C$gene_id <- mapply(strsplit(as.character(m5C$gene_id),","),FUN=function(x){(as.character(x)[1])})
m5C$gene_id <- gsub("Ensembl\\:","",m5C$gene_id)
m5C$m5C_gene_region <- NULL
m5C <- ddply(m5C,"gene_id", numcolwise(sum))
dim(m5C) #4220 genes with m5C sites



### AtoI ###
AtoI <- read.delim("./RMVar_Human_basic_info_A-to-I.txt",sep="\t")
AtoI <- AtoI[!AtoI$gene=="-",]
AtoI <- AtoI[grepl("Ensembl:",AtoI$gene_database_id)==TRUE,]
AtoI <- AtoI[!grepl("intron",AtoI$gene_region)==TRUE,]
AtoI <- AtoI[grepl("Protein coding",AtoI$gene_type)==TRUE,]
AtoI <- data.frame("gene_name"=AtoI$gene,"gene_id"=AtoI$gene_database_id,"gene_region"=AtoI$gene_region,"total_AtoI"=1)
dim(AtoI)

AtoI$CDS_AtoI <- ifelse(grepl("CDS",AtoI$gene_region)==TRUE,1,0)
AtoI$UTR5_AtoI <- ifelse(grepl("5'UTR",AtoI$gene_region)==TRUE,1,0)
AtoI$UTR3_AtoI <- ifelse(grepl("3'UTR",AtoI$gene_region)==TRUE,1,0)
AtoI$gene_id <- mapply(strsplit(as.character(AtoI$gene_id),","),FUN=function(x){(as.character(x)[1])})
AtoI$gene_id <- gsub("Ensembl\\:","",AtoI$gene_id)
AtoI$AtoI_gene_region <- NULL
AtoI <- ddply(AtoI,"gene_id", numcolwise(sum))
dim(AtoI) #1861 genes with AtoI sites




### m1A ###
m1A <- read.delim("./RMVar_Human_basic_info_m1A.txt",sep="\t")
m1A <- m1A[!m1A$gene=="-",]
m1A <- m1A[grepl("Ensembl:",m1A$gene_database_id)==TRUE,]
m1A <- m1A[!grepl("intron",m1A$gene_region)==TRUE,]
m1A <- m1A[grepl("Protein coding",m1A$gene_type)==TRUE,]
m1A <- data.frame("gene_name"=m1A$gene,"gene_id"=m1A$gene_database_id,"gene_region"=m1A$gene_region,"total_m1A"=1)
dim(m1A)

m1A$CDS_m1A <- ifelse(grepl("CDS",m1A$gene_region)==TRUE,1,0)
m1A$UTR5_m1A <- ifelse(grepl("5'UTR",m1A$gene_region)==TRUE,1,0)
m1A$UTR3_m1A <- ifelse(grepl("3'UTR",m1A$gene_region)==TRUE,1,0)
m1A$gene_id <- mapply(strsplit(as.character(m1A$gene_id),","),FUN=function(x){(as.character(x)[1])})
m1A$gene_id <- gsub("Ensembl\\:","",m1A$gene_id)
m1A$m1A_gene_region <- NULL

m1A <- ddply(m1A,"gene_id", numcolwise(sum))
dim(m1A) #7219 genes with m1A sites




### m7G ###
m7G <- read.delim("./RMVar_Human_basic_info_m7G.txt",sep="\t")
m7G <- m7G[!m7G$gene=="-",]
m7G <- m7G[grepl("Ensembl:",m7G$gene_database_id)==TRUE,]
m7G <- m7G[!grepl("intron",m7G$gene_region)==TRUE,]
m7G <- m7G[grepl("Protein coding",m7G$gene_type)==TRUE,]
m7G <- data.frame("gene_name"=m7G$gene,"gene_id"=m7G$gene_database_id,"gene_region"=m7G$gene_region,"total_m7G"=1)
dim(m7G)

m7G$CDS_m7G <- ifelse(grepl("CDS",m7G$gene_region)==TRUE,1,0)
m7G$UTR5_m7G <- ifelse(grepl("5'UTR",m7G$gene_region)==TRUE,1,0)
m7G$UTR3_m7G <- ifelse(grepl("3'UTR",m7G$gene_region)==TRUE,1,0)
m7G$gene_id <- mapply(strsplit(as.character(m7G$gene_id),","),FUN=function(x){(as.character(x)[1])})
m7G$gene_id <- gsub("Ensembl\\:","",m7G$gene_id)
m7G$m7G_gene_region <- NULL

m7G <- ddply(m7G,"gene_id", numcolwise(sum))
dim(m7G) #6281 genes with m7G sites




## conservation ## 
Gene_conservation <- getBM(attributes=c("ensembl_transcript_id", "drerio_homolog_perc_id_r1"), mart = ensembl, filters = "with_ccds", values = TRUE)


## Aggregation of all info into 1 library ##
#UTR3_fasta, CDS_fasta, UTR5_fasta, m6A, m5C, AtoI, m1A, m7G

# importing df #
UTR5_fasta <- read.delim("./5UTR_library_v2_RBP_GC_length.csv",sep=";", dec=",")
CDS_fasta <- read.delim("./CDS_library_v2_RBP_GC_length.csv",sep=";", dec=",")
UTR3_fasta <- read.delim("./3UTR_library_v2_RBP_GC_length.csv",sep=";", dec=",")
miRDB_CD8 <- read.delim("./miRDB_CD8_wrangled_avg_per_tx.csv", sep = ";", dec=",")
Gene_conservation <- Gene_conservation[!duplicated(Gene_conservation$ensembl_transcript_id),]

# formatting #
UTR5_fasta$tx.id <- UTR5_fasta$tx.id_UTR5
CDS_fasta$tx.id <- CDS_fasta$tx.id_CDS
UTR3_fasta$tx.id <- UTR3_fasta$tx.id_UTR3

UTR5_fasta$tx.id_UTR5 <- NULL
CDS_fasta$tx.id_CDS <- NULL
UTR3_fasta$tx.id_UTR3 <- NULL

UTR5_fasta$sequence <- NULL
CDS_fasta$sequence <- NULL
UTR3_fasta$sequence <- NULL

table(duplicated(UTR5_fasta$tx.id))
table(duplicated(CDS_fasta$tx.id))
table(duplicated(UTR3_fasta$tx.id))

# merging df together #
Library_V2 <- NULL
Library_V2 <- merge(UTR5_fasta,CDS_fasta,by="tx.id",all=T)
Library_V2 <- merge(Library_V2,UTR3_fasta,by="tx.id",all=T)
Library_V2 <- merge(Library_V2,IDs_genenames_coding,by.x="tx.id",by.y="ensembl_transcript_id",all.x =T)
Library_V2 <- merge(Library_V2,m6A,by.x="ensembl_gene_id",by.y="gene_id",all.x=T)
Library_V2 <- merge(Library_V2,m5C,by.x="ensembl_gene_id",by.y="gene_id",all.x=T)
Library_V2 <- merge(Library_V2,AtoI,by.x="ensembl_gene_id",by.y="gene_id",all.x=T)
Library_V2 <- merge(Library_V2,m1A,by.x="ensembl_gene_id",by.y="gene_id",all.x=T)
Library_V2 <- merge(Library_V2,m7G,by.x="ensembl_gene_id",by.y="gene_id",all.x=T)
Library_V2 <- merge(Library_V2,Gene_conservation,by.x="tx.id",by.y="ensembl_transcript_id",all.x=T)
Library_V2 <- merge(Library_V2,miRDB_CD8,by.x="tx.id",by.y="ensembl_transcript_id",all.x=T)
# dim(Library_V2)

# Final formatting #
Library_V2 <- subset(Library_V2,Library_V2$transcript_biotype=="protein_coding")
Library_V2$gene_name <- Library_V2$external_gene_name
Library_V2$external_gene_name <- NULL


# Exporting the RNA library ##
write.table(Library_V2,"./RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB.csv",sep = ";", dec = ",",row.names = F)


###__________________________________________________________________________________________###
###______________________________________protein parameters__________________________________###
###__________________________________________________________________________________________###



## Getting Unitprot annotation ##
# Downloaded from Uniprot.org on June 4th 2020
Uniprot_ref <- GeneID_uniprot
Uniprot_ref <- Uniprot_ref[Uniprot_ref$uniprotswissprot>1,]




### PTM annotation ###
# # Downloaded from http://dbptm.mbc.nctu.edu.tw/ on June 4th 2020 # #
# # Only modification found in more than 1000 proteins were acquired # #

## Acetylation ##
Acetylation <- read.delim("./Acetylation.txt",sep = "\t",header = F)
colnames(Acetylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Acetylation$Species <- mapply(strsplit(as.character(Acetylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Acetylation <- subset(Acetylation,Acetylation$Species=="HUMAN")
Acetylation <- data.frame("Entry"=Acetylation$Entry,"Acetylation_count"=1)
Acetylation <- ddply(Acetylation,"Entry", numcolwise(sum))


## Amidation ##
Amidation <- read.delim("./Amidation.txt",sep = "\t",header = F)
colnames(Amidation) <- c("Species","Entry","position","PTM_type","X","Motif")
Amidation$Species <- mapply(strsplit(as.character(Amidation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Amidation <- subset(Amidation,Amidation$Species=="HUMAN")
Amidation <- data.frame("Entry"=Amidation$Entry,"Amidation_count"=1)
Amidation <- ddply(Amidation,"Entry", numcolwise(sum))

## Hydroxylation ##
Hydroxylation <- read.delim("./Hydroxylation.txt",sep = "\t",header = F)
colnames(Hydroxylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Hydroxylation$Species <- mapply(strsplit(as.character(Hydroxylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Hydroxylation <- subset(Hydroxylation,Hydroxylation$Species=="HUMAN")
Hydroxylation <- data.frame("Entry"=Hydroxylation$Entry,"Hydroxylation_count"=1)
Hydroxylation <- ddply(Hydroxylation,"Entry", numcolwise(sum))

## Malonylation ##
Malonylation <- read.delim("./Malonylation.txt",sep = "\t",header = F)
colnames(Malonylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Malonylation$Species <- mapply(strsplit(as.character(Malonylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Malonylation <- subset(Malonylation,Malonylation$Species=="HUMAN")
Malonylation <- data.frame("Entry"=Malonylation$Entry,"Malonylation_count"=1)
Malonylation <- ddply(Malonylation,"Entry", numcolwise(sum))

## Methylation ##
Methylation <- read.delim("./Methylation.txt",sep = "\t",header = F)
colnames(Methylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Methylation$Species <- mapply(strsplit(as.character(Methylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Methylation <- subset(Methylation,Methylation$Species=="HUMAN")
Methylation <- data.frame("Entry"=Methylation$Entry,"Methylation_count"=1)
Methylation <- ddply(Methylation,"Entry", numcolwise(sum))

## N_linked_Glycosylation ##
N_linked_Glycosylation <- read.delim("./N-linkedGlycosylation.txt",sep = "\t",header = F)
colnames(N_linked_Glycosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
N_linked_Glycosylation$Species <- mapply(strsplit(as.character(N_linked_Glycosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
N_linked_Glycosylation <- subset(N_linked_Glycosylation,N_linked_Glycosylation$Species=="HUMAN")
N_linked_Glycosylation <- data.frame("Entry"=N_linked_Glycosylation$Entry,"N_linked_Glycosylation_count"=1)
N_linked_Glycosylation <- ddply(N_linked_Glycosylation,"Entry", numcolwise(sum))

## O_linked_Glycosylation ##
O_linked_Glycosylation <- read.delim("./O-linkedGlycosylation.txt",sep = "\t",header = F)
colnames(O_linked_Glycosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
O_linked_Glycosylation$Species <- mapply(strsplit(as.character(O_linked_Glycosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
O_linked_Glycosylation <- subset(O_linked_Glycosylation,O_linked_Glycosylation$Species=="HUMAN")
O_linked_Glycosylation <- data.frame("Entry"=O_linked_Glycosylation$Entry,"O_linked_Glycosylation_count"=1)
O_linked_Glycosylation <- ddply(O_linked_Glycosylation,"Entry", numcolwise(sum))

## Palmitoylation ##
Palmitoylation <- read.delim("./Palmitoylation.txt",sep = "\t",header = F)
colnames(Palmitoylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Palmitoylation$Species <- mapply(strsplit(as.character(Palmitoylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Palmitoylation <- subset(Palmitoylation,Palmitoylation$Species=="HUMAN")
Palmitoylation <- data.frame("Entry"=Palmitoylation$Entry,"Palmitoylation_count"=1)
Palmitoylation <- ddply(Palmitoylation,"Entry", numcolwise(sum))


## Phosphorylation ##
Phosphorylation <- read.delim("./Phosphorylation.txt",sep = "\t",header = F)
colnames(Phosphorylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Phosphorylation$Species <- mapply(strsplit(as.character(Phosphorylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Phosphorylation <- subset(Phosphorylation,Phosphorylation$Species=="HUMAN")
Phosphorylation <- data.frame("Entry"=Phosphorylation$Entry,"Phosphorylation_count"=1)
Phosphorylation <- ddply(Phosphorylation,"Entry", numcolwise(sum))

## S_nitrosylation ##
S_nitrosylation <- read.delim("./S-nitrosylation.txt",sep = "\t",header = F)
colnames(S_nitrosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
S_nitrosylation$Species <- mapply(strsplit(as.character(S_nitrosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
S_nitrosylation <- subset(S_nitrosylation,S_nitrosylation$Species=="HUMAN")
S_nitrosylation <- data.frame("Entry"=S_nitrosylation$Entry,"S_nitrosylation_count"=1)
S_nitrosylation <- ddply(S_nitrosylation,"Entry", numcolwise(sum))

## Succinylation ##
Succinylation <- read.delim("./Succinylation.txt",sep = "\t",header = F)
colnames(Succinylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Succinylation$Species <- mapply(strsplit(as.character(Succinylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Succinylation <- subset(Succinylation,Succinylation$Species=="HUMAN")
Succinylation <- data.frame("Entry"=Succinylation$Entry,"Succinylation_count"=1)
Succinylation <- ddply(Succinylation,"Entry", numcolwise(sum))

## Sumoylation ##
Sumoylation <- read.delim("./Sumoylation.txt",sep = "\t",header = F)
colnames(Sumoylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Sumoylation$Species <- mapply(strsplit(as.character(Sumoylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Sumoylation <- subset(Sumoylation,Sumoylation$Species=="HUMAN")
Sumoylation <- data.frame("Entry"=Sumoylation$Entry,"Sumoylation_count"=1)
Sumoylation <- ddply(Sumoylation,"Entry", numcolwise(sum))

## Ubiquitination ##
Ubiquitination <- read.delim("./Ubiquitination.txt",sep = "\t",header = F)
colnames(Ubiquitination) <- c("Species","Entry","position","PTM_type","X","Motif")
Ubiquitination$Species <- mapply(strsplit(as.character(Ubiquitination$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Ubiquitination <- subset(Ubiquitination,Ubiquitination$Species=="HUMAN")
Ubiquitination <- data.frame("Entry"=Ubiquitination$Entry,"Ubiquitination_count"=1)
Ubiquitination <- ddply(Ubiquitination,"Entry", numcolwise(sum))


### Combining all annotations and library together ###
Prot_annotation <- merge(Acetylation, Amidation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Hydroxylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Malonylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Methylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, N_linked_Glycosylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, O_linked_Glycosylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Palmitoylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Phosphorylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, S_nitrosylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Succinylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Sumoylation, by="Entry",all=T)
Prot_annotation <- merge(Prot_annotation, Ubiquitination, by="Entry",all=T)

Library_V2_Prot <- merge(Prot_annotation, Uniprot_ref, by.x="Entry",by.y = "uniprotswissprot",all.y=T)

dim(Prot_annotation)
dim(Library_V2_Prot)

Library_V2_Prot[is.na(Library_V2_Prot)]<-0

Library_V2_Prot_per_uniprot_entry <- merge(Library_V2_Prot,Library_V2,by.x="ensembl_transcript_id", by.y="tx.id", all.x=T)

# Dropping unused columns #
#Library_V2_Prot_per_uniprot_entry$ensembl_gene_id.x <- NULL
Library_V2_Prot_per_uniprot_entry$ensembl_gene_id.y <- NULL
Library_V2_Prot_per_uniprot_entry$gene_biotype.x <- NULL
Library_V2_Prot_per_uniprot_entry$Ensembl.transcript <- NULL
Library_V2_Prot_per_uniprot_entry$gene_biotype.y <- NULL
Library_V2_Prot_per_uniprot_entry$external_gene_name <- NULL

## Here I aggregate the features per genes using the mean of all mRNA isoform ##
doMC::registerDoMC(cores = 12) ## Here I use doMC to paralellize the ddply function below
gc()
Library_V2_Prot_per_gene <- plyr::ddply(Library_V2_Prot_per_uniprot_entry,"gene_name", numcolwise(mean),.parallel = T, .progress = T)

write.table(Library_V2_Prot_per_gene,"./Protein_per_gene_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv",sep = "\t",row.names = F)


## Here I aggregate the features per uniprot ID using the mean of all mRNA isoform ##
gc()
Library_V2_Prot_per_uniprot_entry <- plyr::ddply(Library_V2_Prot_per_uniprot_entry,"Entry", numcolwise(mean),.parallel = T, .progress = T)
doMC::registerDoMC(cores = 1) # returning to single core computation

write.table(Library_V2_Prot_per_uniprot_entry,"./Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv",sep = "\t",row.names = F)
# dim(Library_V2_Prot_per_gene)
# dim(Library_V2_Prot_per_uniprot_entry)














