#Example code for TCGA FireBrowser data analysis
rm(list = ls())
#INSTALL PACKAGES Once BY THIS CODE
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
library(stringr) #Run every time
library(limma)
library(edgeR)
library(baySeq)
getwd()
setwd("J:/V/R")

#STAGE-1
a<-read.csv("clinical stage 1234 - cesc.csv")# 16 letters column name
colnames(a)<-toupper(names(a)) # how we have to change col names to upper
colnames(a)<-substr(colnames(a), start = 1, stop = 15)
b<-read.csv("Cancer RSEM.csv") #never use row.names = 1 for merging files
colnames(b)<-substr(colnames(b), start = 1, stop = 16) 
any(is.double(names(a)))
any(is.double(names(b)))
any(which(colnames(b)%in% colnames(a)))
which(colnames(b)%in% colnames(a))
colnames(a) %in% colnames(b) %>% table
c<-merge(a,b,all.x = TRUE,all.y = TRUE)
write.csv(c, "cese RNAseq Stage.csv")
## Change names roman numbers to 1,2,3,4 and note down in number of samples per stage finally uplode
stage<-read.csv("cese RNAseq Stage.csv",header = T, row.names = 1)
names(stage)
stage1<-stage[str_detect(names(stage), "stage.1")]
stage2<-stage[str_detect(names(stage), "stage.2")]
stage3<-stage[str_detect(names(stage), "stage.3")]
stage4<-stage[str_detect(names(stage), "stage.4")]  
stage12<-cbind(stage1, stage2, stage3, stage4) 
index0=floor(dim(stage12)[2]*.2); index0 # Column dimension 
a<-stage12[rowSums(stage1==0)<=index0,]; dim(a) #removing >20% zeros in m0

stage34<-cbind(stage3, stage4) 
index1=floor(dim(stage34)[2]*.2); index1 # Column dimension 
b<-stage34[rowSums(stage34==0)<=index1,]; dim(b) #removing >20% zeros in m0
S<-merge(a, b, by = "row.names", all.x = FALSE,all.y = FALSE); dim(S) # Only matched content will save in data
write.csv(stage12,"cesc RNAseq.csv")
#OR
stage<-read.csv("Cancer RSEM.csv", row.names = 1)
names(stage)
#Stage-1: 45 samples
stage1<-stage[str_detect(names(stage), "stage.1")]
stage1[1,]
colnames(stage1)<-stage1[1,]
stage1<-stage1[-1,]
colnames(stage1) <- paste(colnames(stage1), "S1", sep = "_")
index0=floor(dim(stage1)[2]*.2); index0 # Column dimension 
a<-stage1[rowSums(stage1==0)<=index0,]; dim(a) #removing >20% zeros in m0

#Stage-2: 110 samples
stage2<-stage[str_detect(names(stage), "stage.2")]
stage2[1,]
colnames(stage2)<-stage2[1,]
stage2<-stage2[-1,]
colnames(stage2) <- paste(colnames(stage2), "S2", sep = "_")
index1=floor(dim(stage2)[2]*.2); index1 # Column dimension 
b<-stage2[rowSums(stage2==0)<=index1,]; dim(b) #removing >20% zeros in m0

#Stage-3: 80 samples
stage3<-stage[str_detect(names(stage), "stage.3")]
stage3[1,]
colnames(stage3)<-stage3[1,]
stage3<-stage3[-1,]
colnames(stage3) <- paste(colnames(stage3), "S3", sep = "_")
index2=floor(dim(stage3)[2]*.2); index2 # Column dimension 
c<-stage3[rowSums(stage3==0)<=index2,]; dim(c) #removing >20% zeros in m0

#Stage-4: 39 samples
stage4<-stage[str_detect(names(stage), "stage.4")]  
stage4[1,]
colnames(stage4)<-stage4[1,]
stage4<-stage4[-1,]
colnames(stage4) <- paste(colnames(stage4), "S4", sep = "_")
index3=floor(dim(stage4)[2]*.2); index3 # Column dimension 
d<-stage4[rowSums(stage4==0)<=index3,]; dim(d) #removing >20% zeros in m0
stage12<-cbind(stage1, stage2) 
stage34<-cbind(stage3, stage4) 
stage1234<-cbind(stage12,stage34) 
write.csv(stage1234,"UCEC 20 miR.csv")
#stage12 merge
S12<-merge(a, b, by = "row.names", all.x = FALSE,all.y = FALSE); dim(S12) # Only matched content will save in data
S12[,1] # after merrging we can see row names is in first column 
rownames(S12)<-S12[,1] # we have to replace row names 
S12[,1]
View(S12)
S12<-S12[,-1]
#stage34 merge
S34<-merge(c, d, by = "row.names", all.x = FALSE,all.y = FALSE); dim(S34) # Only matched content will save in data
S34[,1] # after merrging we can see row names is in first column 
rownames(S34)<-S34[,1] # we have to replace row names 
S34[,1]
View(S34)
S34<-S34[,-1]
S1234<-merge(S12, S34, by = "row.names", all.x = FALSE,all.y = FALSE); dim(S1234)
S1234[,1] # after merrging we can see row names is in first column 
rownames(S1234)<-S1234[,1] # we have to replace row names 
S1234[,1]
View(S1234)
S1234<-S1234[,-1]
write.csv(S1234,"UCEC 20 miR.csv")
View(S1234)

#STAGE-3
#DIFFRENTIALLY EXPRESSED ANALYSIS
library(limma)
library(edgeR)
dat1<- read.csv("cesc RNAseq.csv", header = T, row.names = 1) # add First row  Gene name in excel
dim(dat1)
class(dat1)
dat2<-log2(dat1+ 1) # CHANGE DATA TO LOG VALUES 
dat3<-as.matrix(dat1) #CHANGE data frame to matrix form
dim(dat3)
# below is to save design matrix as a separate file and arrange conditions 
write.csv(colnames(dat1),"design matrix.csv") 
design.mat<-read.csv("design matrix.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('S12', 'S34'), "Diff")
sample<-factor(rep(c("S12", "S34"), c(231,66))) #number of T and N samples
design.mat<-model.matrix(~0+sample) 
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = S12 - S34, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat
fit<-lmFit(dat3, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" )  
names(deg1)
class(deg1)
dim(deg1)
write.csv(deg1,"ucec miR DEG.csv")

#STAGE-4
# isolation by using conditions
#condition
UPdeg<-deg1[deg1$logFC>= 1 & deg1$adj.P.Val<0.05,] #UPREGULATED GENES
Downdeg<-deg1[deg1$logFC<= -1 & deg1$adj.P.Val<0.05,] #DOWNREGULATED GENES
DEG_final<-rbind(UPdeg,Downdeg) #5-10
rownames(DEG_final)
library(stringr)
b<-merge(dat1,DEG_final, by = "row.names", all.x = FALSE,all.y = FALSE)
b[,1] # after merging we can see row names is in first column 
rownames(b)<-b[,1] # we have to replace row names 
b[,1]
View(b)
b<-b[,-1]
new_b<-b[,1:176]
bt<-t(new_b)
rownames(bt)
write.csv(rownames(bt),"Class Label.csv") 
# arrange in excel change patient names as row names 
# add new column Class Label Stage1,2 as negative "0" Stage3,4 as Positive "1"
CL<-read.csv("Class Label.csv", row.names = 1)
View(CL)
class(CL)
d<-cbind(bt,CL) #final file for Machine Learning analysis
colnames(d)
class(d)
write.csv(d, "PAAD for SVM.csv")



