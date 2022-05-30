#Example code for TCGA FireBrowser data analysis
# Step-1:
# We have two data files 1.KIRC 20K genes.csv and 2. KIRC DEGs genes for maintain original files safe.
# For first file KIRC 20K genes file copy data and save it as a Cancer_KIRC_20.csv file
# for second file DEGs, use options in .CSV file
#  a) Sort & Filter
#  b) padj Column: Number filter, select Less than or equals to 0.05/0.001/0.0001....etc.
#  c) log2FoldChange: Number filter, select between option greater than or equals to 1 and Less than or equals to -1
#  Here i used padj Less than or equals to 0.0001 and log2FoldChange Absolute |2|, a total of 48 genes
# Copy selected names of genes or Row.names into new file clinical m0 and m1.csv
# finally we make two files i) Cancer_KIRC_20.csv and ii) clinical m0 and m1.csv upload using below code

rm(list = ls())
setwd("J:/V/R")
getwd()
#INSTALL PACKAGES BY THIS CODE
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
library(stringr)
library(limma)
library(edgeR)
getwd()
setwd("J:/V/R")
#Processing m0 and m1 samples
a<-read.csv("Cancer data-2.csv", row.names = 1)
b<-read.csv("Cancer Exp Data.csv", row.names = 1) #use row.names = 1 for merging files
#merging a and b by matching row names 
m0m1<-merge(a,b, by = "row.names", all.x = FALSE,all.y = FALSE); dim(m0m1) # Only matched content will save in data
write.csv(m0m1, "READ.csv") #saving final results
View(m0m1)
# open file KIRC1.csv and remove first column having numbers and save it
# next we needed to add class labels
# We need to open KIRC Phenotype data and select first column "patient IDs" and "pathologic_M"
# Transpose those selected two columns and save these two columns in to new file name "Class Labels.csv"
# Now save and upload two files 
c<-read.csv("READ.csv")
d<-read.csv("Class Labels.csv")
e<-merge(d,c,all.x = TRUE,all.y = TRUE)
# save e as Kirc_CL.csv
# open and arrange Data according to we want in .csv file
# Replace Class Labels M0=0 and M1=1
write.csv(e, "READ_CL.csv")
