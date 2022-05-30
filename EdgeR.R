#https://github.com/sbs87/edger/blob/master/edgeR.md
library(edgeR)
library(baySeq)
getwd()
setwd("J:/V/R")
counts<- read.csv("ucec miR.csv", row.names = 1) # add First row  Gene name in excel
d0 <- DGEList(counts, group = group)
#Pre-processing Calculating Normalization factors
d0 <- calcNormFactors(d0)
dim(d0)
plotMDS(d0, method="bcv", col=as.numeric(d$samples$group))
plotMDS(d0)
design.mat<-read.csv("design matrix.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('S12', 'S34'), "Diff")
group<-factor(rep(c("S12", "S34"), c(291,120))) #number of T and N groups
design.mat<-model.matrix(~0+group) 
colnames(design.mat)<-levels(group)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = S12 - S34, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat

d1 <- estimateCommonDisp(d0, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
dgeTest <- exactTest(d1)
dgeTest
plotBCV(d1)
write.csv(d1, "DEGs edger ucec.csv")
#GLM estimates of dispersion
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d0,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
dgeTest <- exactTest(d2)
dgeTest
plotBCV(d2)


library(edgeR)
cn.color=‘blue’
tr.color=‘brown’
main=‘MDS Plot for Count Data’
colors=c(rep(cn.color,5),rep(tr.color,5))
plotMDS(dge,main=main,labels=colnames(dge$counts),col=colors,las=1)
normalized.counts=cpm(dge)
transposed=t(normalized.counts)
distance=dist(transposed)
clusters=hclust(distance)
plot(clusters)