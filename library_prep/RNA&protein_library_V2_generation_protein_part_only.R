
#title: "ARE_in_hg38"
#author: "Benoit Nicolet"
#date: "2/9/2021"


library(plyr)
library(dplyr)
library(seqinr)
library(stringr)
library(qdapRegex)
library(dplyr)
library(Biostrings)
library(biomaRt)
library(ggplot2)
library(stringdist)
library(fuzzyjoin)
library(coRdon)

# BiocManager::install("coRdon")


setwd("/home/ben/Analysis/RF_human/Library_V2_Aug2021/")



ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#All_filters <- listFilters(ensembl)
#View(All_filters)
#attributes_ens <- listAttributes(ensembl)
#View(attributes_ens)

IDs_genenames_coding <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","transcript_biotype"), mart = ensembl)
IDs_genenames_coding <- IDs_genenames_coding[IDs_genenames_coding$transcript_biotype=="protein_coding",]
dim(IDs_genenames_coding)
#IDs_genenames_coding$transcript_biotype <- NULL

GeneID_uniprot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","uniprotswissprot","transcript_biotype"), mart = ensembl)
GeneID_uniprot <- GeneID_uniprot[GeneID_uniprot$transcript_biotype=="protein_coding",]
GeneID_uniprot$transcript_biotype <- NULL
dim(GeneID_uniprot)




Library_V2 <- read.delim("/home/ben/Analysis/RF_human/Library_V2_Aug2021/RNA_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB.csv",sep = ";", dec = ",")


###__________________________________________________________________________________________###
###______________________________________protein parameters__________________________________###
###__________________________________________________________________________________________###


## Getting Unitprot annotation ##
# Downloaded from Uniprot.org on June 4th 2020
Uniprot_ref <- GeneID_uniprot
Uniprot_ref <- Uniprot_ref[Uniprot_ref$uniprotswissprot>1,]


## Acetylation ##
Acetylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Acetylation.txt",sep = "\t",header = F)
colnames(Acetylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Acetylation$Species <- mapply(strsplit(as.character(Acetylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Acetylation <- subset(Acetylation,Acetylation$Species=="HUMAN")
Acetylation <- data.frame("Entry"=Acetylation$Entry,"Acetylation_count"=1)
Acetylation <- ddply(Acetylation,"Entry", numcolwise(sum))


## Amidation ##
Amidation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Amidation.txt",sep = "\t",header = F)
colnames(Amidation) <- c("Species","Entry","position","PTM_type","X","Motif")
Amidation$Species <- mapply(strsplit(as.character(Amidation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Amidation <- subset(Amidation,Amidation$Species=="HUMAN")
Amidation <- data.frame("Entry"=Amidation$Entry,"Amidation_count"=1)
Amidation <- ddply(Amidation,"Entry", numcolwise(sum))

## Hydroxylation ##
Hydroxylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Hydroxylation.txt",sep = "\t",header = F)
colnames(Hydroxylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Hydroxylation$Species <- mapply(strsplit(as.character(Hydroxylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Hydroxylation <- subset(Hydroxylation,Hydroxylation$Species=="HUMAN")
Hydroxylation <- data.frame("Entry"=Hydroxylation$Entry,"Hydroxylation_count"=1)
Hydroxylation <- ddply(Hydroxylation,"Entry", numcolwise(sum))

## Malonylation ##
Malonylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Malonylation.txt",sep = "\t",header = F)
colnames(Malonylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Malonylation$Species <- mapply(strsplit(as.character(Malonylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Malonylation <- subset(Malonylation,Malonylation$Species=="HUMAN")
Malonylation <- data.frame("Entry"=Malonylation$Entry,"Malonylation_count"=1)
Malonylation <- ddply(Malonylation,"Entry", numcolwise(sum))

## Methylation ##
Methylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Methylation.txt",sep = "\t",header = F)
colnames(Methylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Methylation$Species <- mapply(strsplit(as.character(Methylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Methylation <- subset(Methylation,Methylation$Species=="HUMAN")
Methylation <- data.frame("Entry"=Methylation$Entry,"Methylation_count"=1)
Methylation <- ddply(Methylation,"Entry", numcolwise(sum))

## N_linked_Glycosylation ##
N_linked_Glycosylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/N-linkedGlycosylation.txt",sep = "\t",header = F)
colnames(N_linked_Glycosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
N_linked_Glycosylation$Species <- mapply(strsplit(as.character(N_linked_Glycosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
N_linked_Glycosylation <- subset(N_linked_Glycosylation,N_linked_Glycosylation$Species=="HUMAN")
N_linked_Glycosylation <- data.frame("Entry"=N_linked_Glycosylation$Entry,"N_linked_Glycosylation_count"=1)
N_linked_Glycosylation <- ddply(N_linked_Glycosylation,"Entry", numcolwise(sum))

## O_linked_Glycosylation ##
O_linked_Glycosylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/O-linkedGlycosylation.txt",sep = "\t",header = F)
colnames(O_linked_Glycosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
O_linked_Glycosylation$Species <- mapply(strsplit(as.character(O_linked_Glycosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
O_linked_Glycosylation <- subset(O_linked_Glycosylation,O_linked_Glycosylation$Species=="HUMAN")
O_linked_Glycosylation <- data.frame("Entry"=O_linked_Glycosylation$Entry,"O_linked_Glycosylation_count"=1)
O_linked_Glycosylation <- ddply(O_linked_Glycosylation,"Entry", numcolwise(sum))

## Palmitoylation ##
Palmitoylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Palmitoylation.txt",sep = "\t",header = F)
colnames(Palmitoylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Palmitoylation$Species <- mapply(strsplit(as.character(Palmitoylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Palmitoylation <- subset(Palmitoylation,Palmitoylation$Species=="HUMAN")
Palmitoylation <- data.frame("Entry"=Palmitoylation$Entry,"Palmitoylation_count"=1)
Palmitoylation <- ddply(Palmitoylation,"Entry", numcolwise(sum))


## Phosphorylation ##
Phosphorylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Phosphorylation.txt",sep = "\t",header = F)
colnames(Phosphorylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Phosphorylation$Species <- mapply(strsplit(as.character(Phosphorylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Phosphorylation <- subset(Phosphorylation,Phosphorylation$Species=="HUMAN")
Phosphorylation <- data.frame("Entry"=Phosphorylation$Entry,"Phosphorylation_count"=1)
Phosphorylation <- ddply(Phosphorylation,"Entry", numcolwise(sum))

## S_nitrosylation ##
S_nitrosylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/S-nitrosylation.txt",sep = "\t",header = F)
colnames(S_nitrosylation) <- c("Species","Entry","position","PTM_type","X","Motif")
S_nitrosylation$Species <- mapply(strsplit(as.character(S_nitrosylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
S_nitrosylation <- subset(S_nitrosylation,S_nitrosylation$Species=="HUMAN")
S_nitrosylation <- data.frame("Entry"=S_nitrosylation$Entry,"S_nitrosylation_count"=1)
S_nitrosylation <- ddply(S_nitrosylation,"Entry", numcolwise(sum))

## Succinylation ##
Succinylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Succinylation.txt",sep = "\t",header = F)
colnames(Succinylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Succinylation$Species <- mapply(strsplit(as.character(Succinylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Succinylation <- subset(Succinylation,Succinylation$Species=="HUMAN")
Succinylation <- data.frame("Entry"=Succinylation$Entry,"Succinylation_count"=1)
Succinylation <- ddply(Succinylation,"Entry", numcolwise(sum))

## Sumoylation ##
Sumoylation <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Sumoylation.txt",sep = "\t",header = F)
colnames(Sumoylation) <- c("Species","Entry","position","PTM_type","X","Motif")
Sumoylation$Species <- mapply(strsplit(as.character(Sumoylation$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Sumoylation <- subset(Sumoylation,Sumoylation$Species=="HUMAN")
Sumoylation <- data.frame("Entry"=Sumoylation$Entry,"Sumoylation_count"=1)
Sumoylation <- ddply(Sumoylation,"Entry", numcolwise(sum))

## Ubiquitination ##
Ubiquitination <- read.delim("/home/ben/Analysis/RF_human/sequences/PTMs_04-06-2020/Ubiquitination.txt",sep = "\t",header = F)
colnames(Ubiquitination) <- c("Species","Entry","position","PTM_type","X","Motif")
Ubiquitination$Species <- mapply(strsplit(as.character(Ubiquitination$Species),"\\_"), FUN=function(x){(as.character(x)[2])})
Ubiquitination <- subset(Ubiquitination,Ubiquitination$Species=="HUMAN")
Ubiquitination <- data.frame("Entry"=Ubiquitination$Entry,"Ubiquitination_count"=1)
Ubiquitination <- ddply(Ubiquitination,"Entry", numcolwise(sum))


### Combining stuff together ###
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

#Library_V2_Prot_per_uniprot_entry$ensembl_gene_id.x <- NULL
Library_V2_Prot_per_uniprot_entry$ensembl_gene_id.y <- NULL
Library_V2_Prot_per_uniprot_entry$gene_biotype.x <- NULL
Library_V2_Prot_per_uniprot_entry$Ensembl.transcript <- NULL
Library_V2_Prot_per_uniprot_entry$gene_biotype.y <- NULL
Library_V2_Prot_per_uniprot_entry$external_gene_name <- NULL


# print("by gene")
# 
# doMC::registerDoMC(cores = 12)
# gc()
# Library_V2_Prot_per_gene <- plyr::ddply(Library_V2_Prot_per_uniprot_entry,"gene_name", numcolwise(mean),.parallel = T, .progress = T)
# 
# write.table(Library_V2_Prot_per_gene,"/home/ben/Analysis/RF_human/Library_V2_Aug2021/Protein_per_gene_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv",sep = "\t",row.names = F)


print("by entries")
gc()
Library_V2_Prot_per_uniprot_entry <- plyr::ddply(Library_V2_Prot_per_uniprot_entry,"Entry", numcolwise(mean),.parallel = T, .progress = T)

write.table(Library_V2_Prot_per_uniprot_entry,"/home/ben/Analysis/RF_human/Library_V2_Aug2021/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv",sep = "\t",row.names = F)
#print(dim(Library_V2_Prot_per_gene))
print(dim(Library_V2_Prot_per_uniprot_entry))


