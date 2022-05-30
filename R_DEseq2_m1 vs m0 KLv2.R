#DEGseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

install.packages("pheatmap")
install.packages("tidyverse")
install.packages("stringi")
install.packages("stringr")
library(pheatmap)
library(tidyverse)
library(DESeq2)
library(stringi)
library(stringr)
#STEP-1
#Processing m0 and m1 samples
setwd("/Users/User/Documents/R/xena/KICH")
#[1]download "phenotype" data from xena
#[2]select two columns: patient ID [column A] and 
#the clinical_M or pathologic_M column
#[3]coppy and paste them in transpose with a new Excel file
#[4]name it "clinical_id_metastasis_status.csv"
a<-read.csv("clinical m0 and m1.csv")# 16 letters column name

#colnames(a)<-toupper(names(a)) # how we have to change col names to upper
#[1]download the HTseq-counts data from Xena
#[2]unzip the file and save it TCGA-KICH.htseq_counts.tsv
b<-read.csv("Cancer RSEM.csv") #never use row.names = 1 for merging files
#colnames(b)<-substr(colnames(b), start = 1, stop = 16) 
any(is.double(names(a)))
any(is.double(names(b)))
any(which(colnames(b)%in% colnames(a)))
which(colnames(b)%in% colnames(a))
length(which(colnames(b)%in% colnames(a)))
c<-merge(a,b,all.x = TRUE,all.y = TRUE)
write.csv(c, "TCGA-KICH.htseq_counts-v2.csv") ; 
#dim_c = dim(c); dim_c
#dim_c[1] # total of rows
#c[dim_c[1], 1:5] # make sure that this row consists of M0, M1
#c[1, 1:5] 
#view(c)
#names(c)
#############  read the following 
#[1] need manual edit the file, move bottom line [m0 vs m1 inform] to 1st row
#[2] scan to the right, find the Ensembl_ID column
#[3] cut this column insert it into column A, delete column B [entry num.]
#[4] scan to the right, find where "clinical_M or pathologic_M" column, 
#[4] and delete all the columns [NA]
#[4] onward including the "clinical/pathologic_M" column
#[5] save the file with name "TCGA-KICH.htseq_counts-v3.csv"
cc =read.csv("TCGA-KICH.htseq_counts-v2.csv")
#view(cc)
#which(colnames(cc)=="pathologic_M" ) # it is 45
#cc[1:5, which(colnames(cc)=="pathologic_M" )]
#cc[1:5, 44:46]
#which(cc[1,]=="Ensembl_ID")
#cc[1:4,1]
#cc[1:4,75]
#cc[,1] = cc[,75]
#dim(cc)
#cc[1:4,1]; 
#cc[1:4,2]
#cc[1:4, 44]
#cc[1:4, 45]
#cc[1:4, 46]
#dd = cc[,1:44]
#dim(dd)
#view(dd)
#write.csv(dd, "TCGA-KICH.htseq_counts-v2.csv")
# open "TCGA-KICH.htseq_counts-v3.csv", 1st column is ENSEMBL_ID
#STEP-2
x<-read.csv("TCGA-KICH.htseq_counts-v2.csv", row.names = 1) # column names changed to m0 & m1
y<-x[str_detect(names(x), "M0")] # selection of m0
z<-x[str_detect(names(x), "M1")] #selection of m1
m1m0<-cbind(z,y)

write.csv(m1m0, "TCGA-KICH-rnaseq-m0m1.csv")
view(m1m0)

# after you have the "design_matrix.csv"
#
#[1]need to set up the design matrix, prepare from 'TCGA-KICH-rnaseq-m0m1.csv'
#[2]open it, the first row is somewhat like 'M1.1', 'M1.1' etc., copy 1st/2nd rows
#[3]"design_matrix.csv" compose of 2 columns: 1st column named: patient ID <= from
# 2nd row in 'TCGA-KICH-rnaseq-m0m1.csv'
# 2nd column in "design_matrix.csv" named 'condition'
# depend on cohort M1 first follow by M0 or M0 follow by M1 =>match TCGA-KICH-rnaseq-m0m1.csv
#[4] save the design matrix file in csv format, named 'design_matrix.csv'
# 
#[5] open the "TCGA-KICH-rnaseq-m0m1.csv", delete the 1st row, save it with the same name
#[6] now, "TCGA-KICH-rnaseq-m0m1.csv" has Ensembl_ID in 1st column, TCGA_ID in 1st row
countData <- read.csv("TCGA-KICH-rnaseq-m0m1.csv", row.names=1)
colData <- read.csv("design matrix.csv", row.names = 1)
# converting log2(counts+1) to counts
COU <- (2^countData)-1 #convert log[counts+1] --> counts
write.csv(COU, "counts_KICH_m0vsm1.csv")
COU<-read.csv("counts_KICH_m0vsm1.csv",header=TRUE, row.names=1)
#removing >80% which are zeros
#index0=floor(dim(countData)[2]*.2); index0 # Column dimension 
#a<-countData[rowSums(countData==0)<=index0,]; dim(a)

dim(colData)
dim(countData); dim(COU)
#any(is.na(countData))
#countData<-na.omit(countData)
class(colData)
class(COU)
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData)) # to match names
#countData <- countData[, rownames(colData)]

#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(round(COU) ,colData , design= ~ condition)
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds) # it takes "a few miuntes to finish"
res <- results(dds)
head(res,3); tail(res,3)
class(res)
write.csv(res, "DEGseq2 DERNAs.csv")
res_sig <- subset(res, padj<.05) #Sort by adjusted p-value and display
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
head(res_lfc, 3)
write.csv(res_lfc,"res_lfc_ACC.csv")
plotMA(res)

#Plotting individual genes
#vsd <- vst(dds) # certain plots, we need to normalize our raw count data
#install.packages("pheatmap")
# below commnad has error message
#pheatmap(assay(vsd)[1, ], cluster_rows=TRUE, show_rownames=TRUE,
#         cluster_cols=FALSE, annotation_col=annot_col)

DErna<-read.csv("res_lfc_ACC.csv", row.names = 1)
#Ensembl_2_gene<-read.csv("gencode.gene.info.v22.tsv", row.names = 1, sep="\t")
Ensembl_2_gene<-read.csv("gencode.gene.info.v22.csv", row.names = 1)
DErna_w_name<-merge(Ensembl_2_gene, DErna, by = "row.names", all.x = FALSE,all.y = FALSE); dim(m1m0) # Only matched content will save in data
head(DErna_w_name)
write.csv(DErna_w_name, "DErna_w_name.csv")
# manual edit for log[FC], adj-p<0.001, protein coding genes
signif=read.csv("DErna_w_name.csv", row.names = 1 )
#[1] open DErna_w_name.csv delete 1st column A with entry num.
#[2] save it with the same name DErna_w_name.csv
signif=read.csv("DErna_w_name.csv", row.names = 1 )
data_for_ML = merge(signif, countData , by = "row.names", all.x = FALSE,all.y = FALSE)
#data_for_ML_t = t(data_for_ML)
write.csv(data_for_ML,"data_for_ML.csv")
#[1] delete 1st column
#[2] delete columns without TCGA-xxx-yyy-.... patient ID
#[3] do transpose in Excel, save it as "data_for_ML_t.csv"
#[4] insert M1 and M0 labels from "design_matrix.csv" in the last column
#[5] change M1 -> 1, M0 -> 0 
#[6] save the file as "data_for_ML_final.csv"
#[6] 
