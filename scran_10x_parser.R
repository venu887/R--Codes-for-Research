# This script is designed to parse the 10X scRNA-seq profile by using the Scran package.
#
# Reference: 
#https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html
# https://bioc.ism.ac.jp/packages/3.7/workflows/vignettes/simpleSingleCell/inst/doc/work-4-misc.html  
#Date: Oct 15 2021
# Author: Jeff C.W.Weng

#!/usr/local/bin/R

# Install the required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("edgeR")


# Load the required packages
library(scran)
library(scater)
library(edgeR)

# Read the desired profiles
options(digits = 3)
dge.tn4909281 <- read10X(
  mtx = "GSM4909281_TN-MH0126-matrix.mtx.gz",
  barcodes = "GSM4909281_TN-MH0126-barcodes.tsv.gz",
  gene = "GSE161529_features.tsv.gz",
  DGEList = TRUE)

# Convert gene aliases to official gene symbols
# Annotation file: https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
annot <- alias2SymbolUsingNCBI(dge.tn4909281$genes$Symbol,
  required.columns = c("GeneID", "Symbol"),
  gene.info.file = "Homo_sapiens.gene_info.gz")
dge.tn4909281$genes <- cbind(dge.tn4909281$genes,
  Official = annot$Symbol,
  GeneID = annot$GeneID)

# Quality control 1:
#   (1) Examine the relationship between the library size and the number of expressed genes.
#       Note: Strong linear relationship indicates good sequencing quality.
#   (2) Examine the proportion of reads from mitochondrial genes in each cell.
#       Note: High expression of mitochondrial genes indicate bad sequencing quality of the cell.
mito.genes <- grep("^mt-", dge.tn4909281$genes$Symbol)
percent.mito.genes <- colSums(dge.tn4909281$counts[mito.genes, ]) / dge.tn4909281$samples$lib.size
num.genes <- colSums(dge.tn4909281$counts != 0)
dge.tn4909281$samples <- cbind(dge.tn4909281$samples,
  percent.mito.genes = percent.mito.genes,
  num.genes = num.genes)
par(mfrow = c(1,2))
plot(dge.tn4909281$samples[,c("lib.size","num.genes")],
  pch = 16,
  cex = 0.7)
plot(dge.tn4909281$samples[,c("lib.size","percent.mito.genes")],
  pch = 16,
  cex = 0.7)

# Quality control 2:
#   (1) Filter out genes expressed in fewer than 1% (total cells).
#   (2) Remove genes with no valid official gene symbols.
#   (3) Remove duplicated genes.
desc.sort <- order(rowSums(dge.tn4909281$counts),
  decreasing = TRUE)
dge.tn4909281 <- dge.tn4909281[desc.sort, ]
criterion.1 <- rowSums(dge.tn4909281$counts > 0) >= ncol(dge.tn4909281)*0.01
criterion.2 <- !is.na(dge.tn4909281$genes$Official)
criterion.3 <- !duplicated(dge.tn4909281$genes$Official)
criteria <- criterion.1 & criterion.2 & criterion.3
table(criteria)
dge.tn4909281 <- dge.tn4909281[criterion.1 & criterion.2 & criterion.3, ]
rownames(dge.tn4909281) <- dge.tn4909281$genes$Official

# Exploratory analysis with scran
#   (1) Pre-clustering for performing the cell-specific normalization:
names(dimnames(dge.tn4909281)) <- NULL
single.cell.obj <- SingleCellExperiment(list(counts = dge.tn4909281$counts))
clust <- quickCluster(single.cell.obj)
table(clust)
#   (2) Normalization for adjusting the cell-specific biases:
#       Compute size factors for all cells by using the deconvolution method (Lun, Bach, and Marioni 2016)
#       and then calculate the normalized log-expression values.
#       Note: Strong linear relationship indicated that size factors are well-correlated
#             with the library sizes for all cells.
single.cell.obj <- computeSumFactors(single.cell.obj,
  clusters = clust)
summary(sizeFactors(single.cell.obj))
libSize <- dge.tn4909281$samples$lib.size
par(mfrow = c(1,1))
plot(libSize/1e2,
  sizeFactors(single.cell.obj),
  log = "xy",
  pch = 16,
  cex = 0.7,
  xlab = "Library size (hundreds)",
  ylab = "Size factor")
single.cell.obj <- logNormCounts(single.cell.obj)

#   (3) Examine of the highly variable genes:
#       Decompose the total variance of each gene by fitting a trend to the endogenous variances (Lun, McCarthy, and Marioni 2016).
#       Subtract the fitted value from the total variance for obtaining the highly variable genes.
#       Note: Each dot represents a gene.
#             The blue curve line represents the fitted mean-variance trend.
#             The top highly variable genes appear above the curve.
dec <- modelGeneVar(single.cell.obj)
plot(dec$mean,
  dec$total,
  xlab = "Mean log-expression",
  ylab = "Variance")
curve(metadata(dec)$trend(x),
  col = "blue",
  add=TRUE)

#   (4) Extraction of top highly variable genes for use in downstream analysis:
#       Get the top 1000 genes.
top.hvgs <- getTopHVGs(dec,
  n = 1000)
single.cell.obj <- fixedPCA(single.cell.obj,
  subset.row = top.hvgs)
reducedDimNames(single.cell.obj)

#   (5) Principle component analysis for retaining the crucial cell clusters:
#       Choose the number of PCs that is not less than the number of clusters.
output <- getClusteredPCs(reducedDim(single.cell.obj))
npcs <- metadata(output)$chosen
reducedDim(single.cell.obj, "PCAsub") <- reducedDim(single.cell.obj, "PCA")[, 1:npcs, drop = FALSE]

#   (6) Graph-based clustering:
#       Cluster the shared nearest neighbors (Xu and Su 2015) and
#       visualize this clustering on a t-SNE plot.
g <- buildSNNGraph(single.cell.obj,
  use.dimred = "PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership
colLabels(single.cell.obj) <- factor(cluster)
table(colLabels(single.cell.obj))
single.cell.obj <- runTSNE(single.cell.obj,
  dimred = "PCAsub")
plotTSNE(single.cell.obj,
  colour_by = "label",
  text_by = "label")

#   (7) Identification of marker genes:
#       Select the top 10 markers of all the clusters.
#       Note: A positive minimum AUC means that the gene is upregulated against all other clusters.
markers <- scoreMarkers(single.cell.obj,
  full.stats = TRUE)
top10 <- lapply(markers, "[", 1:10, )
GSet <- as.vector(sapply(top10, "rownames"))
GSet <- GSet[!duplicated(GSet)]
#       Sort the mean AUC in decreasing order for identifying the upregulated in the cluster 1 compared to other clusters.
markers[[1]][order(markers[[1]]$mean.AUC, decreasing=TRUE),1:19]

#   (8) Detection of the correlated genes
#       Calculate the correlation with the top 1000 genes and another 
cor.pairs <- correlatePairs(single.cell.obj,
  subset.row = top.hvgs)
#cor.pairs <- correlatePairs(single.cell.obj.normal, subset.row = top.hvgs)


#1. Introduction
#https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")
library(scRNAseq)
#2. Setting up the data
sce <- GrunPancreasData()
sce
library(scuttle)
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats, percent_subsets="altexps_ERCC_percent")
sce <- sce[,!qcfilter$discard]
summary(qcfilter$discard)
library(scran)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
summary(sizeFactors(sce))
sce <- logNormCounts(sce)
#3. Variance modelling
dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
dec2 <- modelGeneVarWithSpikes(sce, 'ERCC')
plot(dec2$mean, dec2$total, xlab="Mean log-expression", ylab="Variance")
points(metadata(dec2)$mean, metadata(dec2)$var, col="red")
curve(metadata(dec2)$trend(x), col="blue", add=TRUE)
# Turned off weighting to avoid overfitting for each donor.
dec3 <- modelGeneVar(sce, block=sce$donor, density.weights=FALSE)
per.block <- dec3$per.block
par(mfrow=c(3, 2))
for (i in seq_along(per.block)) {
  decX <- per.block[[i]]
  plot(decX$mean, decX$total, xlab="Mean log-expression", 
       ylab="Variance", main=names(per.block)[i])
  curve(metadata(decX)$trend(x), col="blue", add=TRUE)
}
# Get the top 10% of genes.
top.hvgs <- getTopHVGs(dec, prop=0.1)

# Get the top 2000 genes.
top.hvgs2 <- getTopHVGs(dec, n=2000)

# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(dec, var.threshold=0)

# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)

sce <- fixedPCA(sce, subset.row=top.hvgs)
reducedDimNames(sce)
## [1] "PCA"
#4.Automated PC choice Prinicipal Component
sced <- denoisePCA(sce, dec2, subset.row=getTopHVGs(dec2, prop=0.1))
ncol(reducedDim(sced, "PCA"))
output <- getClusteredPCs(reducedDim(sce))
npcs <- metadata(output)$chosen
reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]
npcs
#5.Graph-based clustering
# In this case, using the PCs that we chose from getClusteredPCs().
g <- buildSNNGraph(sce, use.dimred="PCAsub") #shared nearest neighbors method
cluster <- igraph::cluster_walktrap(g)$membership

# Assigning to the 'colLabels' of the 'sce'.
colLabels(sce) <- factor(cluster)
table(colLabels(sce))
library(scater)
sce <- runTSNE(sce, dimred="PCAsub")
plotTSNE(sce, colour_by="label", text_by="label") #t -SNE plot

library(bluster)
ratio <- pairwiseModularity(g, cluster, as.ratio=TRUE)

library(pheatmap)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))
#A more general diagnostic involves bootstrapping to determine the stability of the partitions between clusters
ass.prob <- bootstrapStability(sce, FUN=function(x) {
  g <- buildSNNGraph(x, use.dimred="PCAsub")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce$cluster)

pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=colorRampPalette(c("white", "blue"))(100))

subout <- quickSubCluster(sce, groups=colLabels(sce))
table(metadata(subout)$subcluster) # formatted as '<parent>.<subcluster>'

#6.Identifying marker genes
#The scoreMarkers() wrapper function will perform differential expression comparisons between pairs of clusters to identify potential marker genes
# Uses clustering information from 'colLabels(sce)' by default:
markers <- scoreMarkers(sce)
markers
colnames(markers[[1]])
# Just showing the first few columns for brevity.
markers[[1]][order(markers[[1]]$mean.AUC, decreasing=TRUE),1:4]
markers <- scoreMarkers(sce, full.stats=TRUE)
markers[[1]]$full.logFC.cohen
#7.Detecting correlated genes
# Using the first 200 HVs, which are the most interesting anyway.
of.interest <- top.hvgs[1:200]
cor.pairs <- correlatePairs(sce, subset.row=of.interest)
cor.pairs
cor.pairs2 <- correlatePairs(sce, subset.row=of.interest, block=sce$donor)
cor.genes <- correlateGenes(cor.pairs)
cor.genes

#8. Converting to other formats
y <- convertTo(sce, type="edgeR")
sessionInfo()
