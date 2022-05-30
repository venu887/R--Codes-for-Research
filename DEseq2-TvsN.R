#STEP-1
#DEGseq2
setwd("J:/V/R")
rm(list = ls())
library(pheatmap)
library(tidyverse)
library(DESeq2)
library(stringr)
# isolation of T and N samples
df1<-read.csv("Cancer RSEM.csv", row.names = 1)
y<-df1[str_detect(colnames(df1), "01A|02A|03A|04A|05A|06A|07A|08A|09A
              |01B|02B|03B|04B|05B|06B|07B|08B|09B
              |01C|02C|03C|04C|05C|06C|07C|08C|09C
              |01D|02D|03D|04D|05D|06D|07D|08D|09D")]

z<-df1[str_detect(names(df1), "11A")]
#merging y and z by matching row names 
TuNo<-merge(y, z, by = "row.names", all.x = FALSE,all.y = FALSE)# Only matched content will save in data
# arrange everything in excel and upload it for further analysis
# i will arrange data by using R
TuNo[1,] #all 1 st columns
TuNo[,1] #all 1 st row
rownames(TuNo)<-TuNo[,1]
View(TuNo)
TuNo<-TuNo[,-1] # Removing first column with gene names 
write.csv(TuNo,"TN 20.csv") #MUST SAVE FILE

#STEP-2
countData <- read.csv("TN 20.csv", row.names=1)
colData <- read.csv("design matrix.csv", row.names = 1)
# converting log2(counts+1) to counts
COU <- (2^countData)-1
write.csv(COU, "ESCA counts.csv")
COU2<-read.csv("ESCA counts.csv", row.names=1)
#removing >80% which are zeros
index0=floor(dim(countData)[2]*.2); index0 # Column dimension 
a<-countData[rowSums(countData==0)<=index0,]; dim(a)

dim(colData)
dim(countData)
any(is.na(COU2))
#countData<-na.omit(countData)
class(colData)
class(countData)
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(COU2)) # to match names
#countData <- countData[, rownames(colData)]

#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(round(COU2), colData, design= ~ condition)
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
write.csv(res, "DEGs DEGseq2 ESCA.csv")
res_sig <- subset(res, padj<.05) #Sort by adjusted p-value and display
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
head(res_lfc)
plotMA(res)
plotCounts(dds, gene=which.min(res$padj), intgroup="Stages.pvm") #Plotting individual genes
vsd <- vst(dds) # certain plots, we need to normalize our raw count data
install.packages("pheatmap")

pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_col)
