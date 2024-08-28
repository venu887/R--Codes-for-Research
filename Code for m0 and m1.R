#Example code for TCGA FireBrowser data analysis
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
rm(b)

#Processing m0 and m1 samples
a<-read.csv("Class Labels.csv", row.names = 1)# 16 letters column name
colnames(a)<-toupper(names(a)) # how we have to change col names to upper
colnames(a)<-substr(colnames(a), start = 1, stop = 16)
b<-read.csv("Cancer Exp Data.csv", row.names = 1) #never use row.names = 1 for merging files
colnames(b)<-substr(colnames(b), start = 1, stop = 16) 
any(is.double(names(a)))
any(is.double(names(b)))
any(which(colnames(b)%in% colnames(a)))
which(colnames(b)%in% colnames(a))
common_cols <- intersect(colnames(a), colnames(b))
colnames(b) %in% colnames(a) %>% table
c<-merge(a,b,all.x = TRUE,all.y = TRUE)
#d<-c[1:20504,1:900] #Extracting large portion to small portion. 
write.csv(c, "Cancer RSEM.2.csv")
# Arrange data and upload .csv file and upload in R
x<-read.csv("Cancer RSEM.2.csv", row.names = 1) # column names changed to m0 & m1
y<-x[str_detect(names(x), "M0")] # selection of m0
write.csv(y, "KIRP_ICGC.csv")
colnames(d1)<-substr(colnames(d1), start = 1, stop = 12)
any(is.double(names(d1)))


z<-x[str_detect(names(x), "M1")]#selection of m1
d2<-b[str_detect(colnames(b), "11A")]
write.csv(d1, "COAD m0 m1.1.csv")
colnames(d2)<-substr(colnames(d2), start = 1, stop = 12)
any(is.double(names(d2)))

any(which(colnames(d1)%in% colnames(d2)))
which(colnames(d1)%in% colnames(d2))
colnames(z) %in% colnames(y) %>% table
m0m1<-cbind(z,y)
write.csv(m0m1, "Cancer RSEM.csv")
#rearrange in excel remove m0m1 and arrange column names as TCGA code and note down m0,m1 number of samples
#STEP-2 removing >20% zeros and rearrange data.
# verifying number of zero in a row and remove samples have more than 20%
# Ex:  (284/100)*20
# if got 2.83*20 like that then take it as 57 samples then extract.
# if total zeros in row we will select 20% by using below code.
setwd("J:/V/R"); list.files()
library(stringr)
library(proxyC) #install this package. 
x<-read.csv("Cancer RSEM.csv", row.names = 1) # column names changed to m0 & m1
#(196/100)*20  we need to know how many samples for 20% for both m0 and m1
y<-x[str_detect(names(x), "m0")]; dim(y) # selection of m0
y[1,]
colnames(y)<-y[1,]#we are changing colnames with TCGA ID along with m0
colnames(y) <- paste(colnames(y), "m0", sep = "_")
View(y)
y<-y[-1,]
View(y)

index0=floor(dim(m0)[2]*.2); index0 # Column dimension 
a<-m0[rowSums(m0==0)<=index0,]; dim(m0) #removing >20% zeros in m0
m0<-b[str_detect(colnames(b), "01A|02A|03A|04A|05A|06A|07A|08A|09A
              |01B|02B|03B|04B|05B|06B|07B|08B|09B
              |01C|02C|03C|04C|05C|06C|07C|08C|09C
              |01D|02D|03D|04D|05D|06D|07D|08D|09D"),]

z<-x[str_detect(names(x), "m1")]; dim(z)#selection of m1
z[1,]
colnames(z)<-z[1,] #we are changing colnames with TCGA ID along with m1
colnames(z) <- paste(colnames(z), "m1", sep = "_")
View(z)
z<-z[-1,]

index1=floor(dim(z)[2]*.2); index1 # Column dimension of z
b<-z[rowSums(z==0)<=index1,]; dim(b) #removing >20% zeros in m1
m1<-z[str_detect(colnames(z), "01A|02A|03A|04A|05A|06A|07A|08A|09A
              |01B|02B|03B|04B|05B|06B|07B|08B|09B
              |01C|02C|03C|04C|05C|06C|07C|08C|09C
              |01D|02D|03D|04D|05D|06D|07D|08D|09D")]

#merging a and b by matching row names 
m0m1<-merge(a,b, by = "row.names", all.x = FALSE,all.y = FALSE); dim(m0m1) # Only matched content will save in data
View(m0m1)
any(is.na(m0m1))
m0m1[,1] # after merrging we can see row names is in first column 
rownames(m0m1)<-m0m1[,1] # we have to replace row names 
m0m1[,1]
View(m0m1)
m0m1<-m0m1[,-1] # we are going to remove extra row names 
write.csv(m0m1, "CCLE 3 cell lines.csv") #saving final results
View(m0m1)


rm(list = ls())
#Generating raw counts from log transformed values 
CO<-read.csv("THCA-20% zeros UCSC m1m0.csv", row.names = 1) # column names changed to m0 & m1
LUAD_COU <- (2^CO)-1
LUAD_COU<-round(LUAD_COU)
write.csv(LUAD_COU, "THCA counts.csv")
#DIFFRENTIALLY EXPRESSED ANALYSIS
library(limma)
library(edgeR)
dat1<-read.csv("Cancer RSEM.csv", row.names = 1)
#d1<-dat1[str_detect(colnames(dat1), "01A")]
any(is.na(dat1))
dim(dat1)
class(dat1)
dat2<-log2(dat1+ 1) # CHANGE DATA TO LOG VALUES 
dat3<-as.matrix(dat1) #CHANGE data frame to matrix form
dim(dat3)
# below is to save design matrix as a separate file and arrange conditions 
write.csv(colnames(dat3),"design matrix m1 m0.csv") 
design.mat<-read.csv("design matrix m1 m0.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('m1', 'm0'), "Diff")
sample<-factor(rep(c("m1", "m0"), c(14,25))) #number of m1 and m0 samples
design.mat<-model.matrix(~0+sample) 
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = m1 - m0, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat
fit<-lmFit(dat3, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" ) # top genes 
names(deg1)
class(deg1)
dim(deg1)
write.csv(deg1,"GEO.D_DEGs.csv")

# isolation by using conditions
UPdeg<-deg1[deg1$logFC>= 0.5 & deg1$adj.P.Val<0.05,] #UPREGULATED GENES
Downdeg<-deg1[deg1$logFC<= -0.5 & deg1$adj.P.Val<0.05,] #DOWNREGULATED GENES
# Data Arrangement like IRIS data for Linear model regression
# upload data that you used for dat1
TotalDEG<-rbind(UPdeg,Downdeg) # must less than 10
View(TotalDEG)
library(stringr)
Cohort_L<-merge(dat1,TotalDEG, by = "row.names", all.x = FALSE,all.y = FALSE)
Cohort_L[,1] # after merging we can see row names is in first column 
rownames(Cohort_L)<-Cohort_L[,1] # we have to replace row names 
Cohort_L[,1]
Cohort_L<-Cohort_L[,-1]
new_C<-Cohort_L[,1:77]# select only RSEM values
Ct<-t(new_C) #transpose data
rownames(Ct)
write.csv(rownames(Ct), "Class Label.csv")
# arrange in excel change patient names as row names 
# add new column Class Label m0 as negative "0" m1 as Positive "1"
CL<-read.csv("Class Label.csv", row.names = 1)
View(CL)
class(CL)
d<-cbind(Ct,CL) #final file for Machine Learning analysis
colnames(d)
class(d)
write.csv(d, "ACC for SVM mRNA.csv") 

indexes<-sample(nrow(d),0.8*nrow(d),replace = F)
indexes
sort(indexes)
train<-d[indexes,]; dim(train)
class(train)
new_train<-as.data.frame(train)
test<-d[-indexes,]; dim(test)
class(test)
test<-as.data.frame(test)
model=glm(c~.,family=binomial,data=new_train)#AIC: 10
summary(model)
model_pred=zapsmall( predict(model, newdata = new_train, type = "response") )
View(model_pred)
model_pred
test_pred = predict(model, newdata = test, type = "response")
test_pred=round(test_pred,1); test_pred
View(test_pred)
length(test_pred)


#miR
Acc_miR<-read.csv("ACC miR.csv", row.names = 1)
Acc_miR[1,]
acc<-Acc_miR[str_detect(Acc_miR[1,], "read_count")]
write.csv(acc, "acc read_count.csv")
a<-read.csv("acc read_count.csv")
colnames(a)<-substr(colnames(a), start = 1, stop = 16) 

b<-read.csv("clinical m0 and m1.csv")# 16 letters column name
colnames(b)<-toupper(names(b)) # how we have to change col names to upper
any(is.double(names(a)))
any(is.double(names(b)))
any(which(colnames(b)%in% colnames(a)))
which(colnames(b)%in% colnames(a))
colnames(b) %in% colnames(a) %>% table
rm(d)
c<-merge(b,a,all.x = TRUE,all.y = TRUE)
write.csv(c, "ACC miR merged.csv")
library(stringr)
library(proxyC) #install this package. 
x<-read.csv("ACC miR merged.csv", row.names = 1) # column names changed to m0 & m1
#(196/100)*20  we need to know how many samples for 20% for both m0 and m1
y<-x[str_detect(names(x), "m0")]; dim(y) # selection of m0
y[1,]
colnames(y)<-y[1,]#we are changing colnames with TCGA ID along with m0
colnames(y) <- paste(colnames(y), "m0", sep = "_")
View(y)
y<-y[-1,]
View(y)
index0=floor(dim(y)[2]*.2); index0 # Column dimension 
a<-y[rowSums(y==0)<=index0,]; dim(a) #removing >20% zeros in m0

z<-x[str_detect(names(x), "m1")]; dim(z)#selection of m1
z[1,]
colnames(z)<-z[1,] #we are changing colnames with TCGA ID along with m1
colnames(z) <- paste(colnames(z), "m1", sep = "_")
View(z)
z<-z[-1,]
index1=floor(dim(z)[2]*.2); index1 # Column dimension of z
b<-z[rowSums(z==0)<=index1,]; dim(b) #removing >20% zeros in m1
#merging a and b by matching row names 
m0m1<-merge(a, b, by = "row.names", all.x = FALSE,all.y = FALSE); dim(m0m1) # Only matched content will save in data
View(m0m1)
any(is.na(m0m1))
m0m1[,1] # after merrging we can see row names is in first column 
rownames(m0m1)<-m0m1[,1] # we have to replace row names 
m0m1[,1]
View(m0m1)
m0m1<-m0m1[,-1] # we are going to remove extra row names 
write.csv(m0m1, "ACC miR20.csv") #saving final results
View(m0m1)


#DIFFRENTIALLY EXPRESSED ANALYSIS
library(limma)
library(edgeR)
dat1<-read.csv("THCA all.csv", row.names = 1)
any(is.na(dat1))
dim(dat1)
class(dat1)
dat2<-log2(dat1+ 1) # CHANGE DATA TO LOG VALUES 
dat3<-as.matrix(dat2) #CHANGE data frame to matrix form
dim(dat3)
# below is to save design matrix as a separate file and arrange conditions 
write.csv(colnames(dat1),"design matrix m0 m1.csv") 
design.mat<-read.csv("design matrix.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('m1', 'm0'), "Diff")
sample<-factor(rep(c("m1", "m0"), c(9,279))) #number of m0 and m1 samples
design.mat<-model.matrix(~0+sample) 
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = m1 - m0, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat
fit<-lmFit(dat3, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" ) # top genes 
names(deg1)
class(deg1)
dim(deg1)
write.csv(deg1,"m1m0 THCA DEG.csv")

# isolation by using conditions
UPdeg<-deg1[deg1$logFC>= 1 & deg1$adj.P.Val<0.05,] #UPREGULATED GENES
Downdeg<-deg1[deg1$logFC<= -1 & deg1$adj.P.Val<0.05,] #DOWNREGULATED GENES

#DEGseq2
countData <- read.csv("THCA m0m1.csv",header=TRUE, row.names=1)
colData <- read.csv("design matrix.csv", row.names = 1)
COU <- (2^countData)-1
write.csv(COU, "THCA counts.csv")
COU<-read.csv("ACC counts.csv",header=TRUE, row.names=1)
index0=floor(dim(CESE_COU)[2]*.2); index0 # Column dimension 
a<-BRCA_COU[rowSums(CESE_COU==0)<=index0,]; dim(a)
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
library(tidyverse)
library(DESeq2)
any(is.integer(countData))
any(is.numeric(countData))
?`DESeq2-package`

#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(round(COU) ,colData , design = ~ Condition)
# filter genes where there are less than 5 samples with normalised counts greater or equal to 10
keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 5
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
write.csv(res, "DEGs DEseq2-20% Zeros removal THCA.csv")
res_sig <- subset(res, padj<.05) #Sort by adjusted p-value and display
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
head(res_lfc)
plotMA(res)
plotCounts(dds, gene=which.min(res$padj), intgroup="Stages.pvm") #Plotting individual genes
vsd <- vst(dds) # certain plots, we need to normalize our raw count data
install.packages("pheatmap")
library(pheatmap)
pheatmap(assay(vsd), cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_col)


