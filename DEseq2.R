#DEGseq2
rm(list = ls())
library(pheatmap)
library(tidyverse)
library(DESeq2)
setwd("J:/V/R")
countData <- read.csv("Cancer RSEM.csv", row.names=1)
write.csv(colnames(countData),"design matrix.csv")
colData <- read.csv("design matrix.csv", row.names = 1)
# converting log2(counts+1) to counts
COU <- (2^countData)-1
write.csv(COU, "Cancer RSEM.3.csv")

COU<-read.csv("Cancer RSEM.3.csv",header=TRUE, row.names=1)
#removing >80% which are zeros
index0=floor(dim(countData)[2]*.2); index0 # Column dimension 
a<-countData[rowSums(countData==0)<=index0,]; dim(a)

str(countData)
rownames(countData)
rownames(colData)
colnames(colData)
colnames(countData)
dim(colData)
dim(countData)

any(is.na(countData))
countData<-na.omit(countData)
class(colData)
class(countData)

#Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData)) # to match names
countData <- countData[, rownames(colData)]
View(countData)
all(rownames(colData) == colnames(countData))

dups <- which(duplicated(rownames(countData)))
length(dups)
unique(countData[dups])
anyDuplicated(rownames(countData))

#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(round(COU) ,colData , design= ~ condition)
#Run the default analysis for DESeq2 and generate results table
#keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
summary(res)
topT <- as.data.frame(res)
head(topT)
names(topT)
write.csv(res, "Cancer RSEM.2.csv")
#Merging Ensymbal data with gene names 
a<-read.csv("Cancer RSEM.2.csv",row.names = 1)
b<-read.csv("Ensymbal.csv",row.names = 1) 
DEGs<-merge(a, b, by = "row.names", all.x = FALSE,all.y = FALSE)
write.csv(DEGs, "Cancer RSEM.2.csv")

#write.csv(topT, "DEGs DEGseq2 CESC mRNA 20.csv")
res_sig <- subset(res, padj<.05) #Sort by adjusted p-value and display
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
# Extract DEGs from Count data
# it has 6 additinal columns from res_lfc 
Expression_of_expression_DEGs<-merge(COU,res_lfc, by = "row.names", all.x = FALSE,all.y = FALSE)


head(res_lfc)
plotMA(res)
resultsNames(dds)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene=which.min(res_sig$padj), intgroup="Condition") #Plotting individual genes
vsd <- vst(dds) # certain plots, we need to normalize our raw count data
#install.packages("pheatmap")

pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_col)
#Volcano plots
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="cesc DEmiR's early vs late", cex=1.0,col="light blue", xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~adj.P~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="purple", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="brown", lty=4, lwd=2.0)
abline(v=1, col="brown", lty=4, lwd=2.0)
alpha <- 0.05 # Threshold on the adjusted p-value
abline(h=-log10(alpha), col="brown")
gn.selected <- abs(res$log2FoldChange) > 1 & topT$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.5)
#if (!requireNamespace('BiocManager', quietly = TRUE))
  #install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(topT,
                lab = rownames(topT),
                x = 'log2FoldChange',
                y = 'pvalue',
                pointSize = 3.0,
                title = 'ov DEmiR early vs late',
                pCutoff = 0.0500,
                FCcutoff = 1,
                shape = c(1, 4, 23, 25),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 2.0,)
## plot using ggplot2
ggplot(topT) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype), position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("miR") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))

# We select gene names based on FDR (1%)
epsilon <- 1 # pseudo-count to avoid problems with log(0)
gene.kept <- rownames(topT)[topT$padj <= alpha & !is.na(topT$padj)]

# We retrieve the normalized counts for gene of interest
countData.kept <- log2(countData + epsilon)[gene.kept, ]
dim(countData.kept)
library("gplots")
heatmap.2(as.matrix(countData.kept), 
          scale="row", 
          hclust=function(x) hclust(x,method="average"), 
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          trace="none", 
          density="none", 
          labRow="",
          cexCol=0.7)
