#Differential Expression with Limma-Voom
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
setwd("J:/V/R")
rm(list = ls())
library(edgeR)
library(limma)
counts <- read.csv("cancer RSEM.csv", row.names = 1)
head(counts)
d0 <- DGEList(counts)
#Pre-processing Calculating Normalization factors
d0 <- calcNormFactors(d0)

#or
dge <- calcNormFactors(d0, method = "TMM")  # normalization
#Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
design.mat<-read.csv("design matrix.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('NS', 'SM'), "Diff")
sample<-factor(rep(c('NS', 'SM'), c(66,81))) #number of T and N samples
design.mat<-model.matrix(~0+sample) 
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = NS - SM, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat
#VOOM plots
y <- voom(d, design.mat, plot = T) # CPM 
tmp <- voom(d0, design.mat, plot = T)
#lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(y, design.mat)
head(coef(fit))
tmp <- contrasts.fit(fit, contrast.mat)
tmp <- eBayes(tmp)

# genes are most differentially expressed
tempOutput = topTable(tmp, n=Inf, adjust.method = 'BH',coef=1)
length(which(tempOutput$adj.P.Val < 0.05))
# genes are most differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
#logFC: log2 fold change of S12/S34
#AveExpr: Average expression across all samples, in log2 CPM
#t: logFC divided by its standard error
#P.Value: Raw p-value (based on t) from test that logFC differs from 0
#adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
#B: log-odds that gene is DE (arguably less useful than the other columns)

#How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))

# This is correct, as we can see from the most upregulated genes.
ix <- match(row.names(top.table)[1], row.names(dat1) )
top1 <- ( as.vector( t( dat1[ix,]) ) )
names( top1 ) = colnames(dat1)
barplot( top1 , las=3)

write.csv(top.table,"PAAD DEG limma-voom.csv")
