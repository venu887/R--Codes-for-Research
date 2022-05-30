#Example code for TCGA FireBrowser data analysis
rm(list = ls())
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

#STEP-1
# isolation of T and N samples
df1<-read.csv("Cancer RSEM.csv", row.names = 1)
y<-df1[str_detect(colnames(df1), "01A|02A|03A|04A|05A|06A|07A|08A|09A
              |01B|02B|03B|04B|05B|06B|07B|08B|09B
              |01C|02C|03C|04C|05C|06C|07C|08C|09C
              |01D|02D|03D|04D|05D|06D|07D|08D|09D")]
index0=floor(dim(y)[2]*.2); index0 # Column dimension 
a<-y[rowSums(y==0)<=index0,]; dim(a) #removing >20% zeros 

z<-df1[str_detect(names(df1), "11A")]
index1=floor(dim(z)[2]*.2); index1 # Column dimension of z
b<-z[rowSums(z==0)<=index1,]; dim(b) #removing >20% zeros 
#merging a and b by matching row names 
TuNo<-merge(a, b, by = "row.names", all.x = FALSE,all.y = FALSE)# Only matched content will save in data
# arrange everything in excel and upload it for further analysis
# i will arrange data by using R
TuNo[1,] #all 1 st columns
TuNo[,1] #all 1 st row
rownames(TuNo)<-TuNo[,1]
View(TuNo)
TuNo<-TuNo[,-1] # Removing extra row names 
write.csv(TuNo,"TN 20.csv") #MUST SAVE FILE
df1234<-cbind(y,z) # row names plus tumour plus normal
write.csv(df1234,"TN all.csv")


#DIFFRENTIALLY EXPRESSED ANALYSIS
library(limma)
library(edgeR)
dat1<-read.csv("E.csv", header = T, row.names = 1)
dim(dat1)
class(dat1)
dat2<-log2(dat1+ 1) # CHANGE DATA TO LOG VALUES 
dat3<-as.matrix(dat2) #CHANGE data frame to matrix form
dim(dat3)
# below is to save design matrix as a separate file and arrange conditions 
write.csv(colnames(dat1),"design matrix TN.csv") # Arrange and add columns T=1 or 0, N=1 or 0
design.mat<-read.csv("design matrix TN.csv") # design of data
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('T', 'N'), "Diff")
sample<-factor(rep(c("T", "N"), c(185,11))) #number of T and N samples
design.mat<-model.matrix(~0+sample) 
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = T - N, levels = design.mat) #jouning of two datas by contrast conditions
contrast.mat
fit<-lmFit(dat3, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" ) # top genes 
dim(deg1)
names(deg1)
class(deg1)
str(deg1)
write.csv(deg1,"TN DEG all.csv")

# Data Arrangement like IRIS data for Linear model regression
# upload data that you used for dat1
UPdeg<-deg1[deg1$logFC>= 6 & deg1$adj.P.Val<0.05,] #UPREGULATED GENES
Downdeg<-deg1[deg1$logFC<= -6 & deg1$adj.P.Val<0.05,] #DOWNREGULATED GENES
TotalDEG<-rbind(UPdeg,Downdeg) # must less than 10
View(TotalDEG)
library(stringr)
Cohort_L<-merge(dat1,TotalDEG, by = "row.names", all.x = FALSE,all.y = FALSE)
Cohort_L[,1] # after merging we can see row names is in first column 
rownames(Cohort_L)<-Cohort_L[,1] # we have to replace row names 
Cohort_L[,1]
Cohort_L<-Cohort_L[,-1]
new_C<-Cohort_L[,1:196]# select only RSEM values
Ct<-t(new_C) #transpose data
# save data and arrange for condition
rownames(Ct)
write.csv(rownames(Ct),"Class Label.csv") 
# arrange in excel change patient names as row names 
# add new column Class Label T as Positive "1" N as negative "0"
CL<-read.csv("Class Label.csv", row.names = 1)
View(CL)
class(CL)
str(CL)
d<-cbind(Ct,CL) #final file for Machine Learning analysis
colnames(d)
class(d)
View(d)

set.seed(2323)
indexes<-sample(nrow(d),0.8*nrow(d),replace = F)
indexes
sort(indexes)
train<-d[indexes,]; dim(train)
class(train)
new_train<-as.data.frame(train)
test<-d[-indexes,]; dim(test)
class(test)
#test<-as.data.frame(test) #if it is already dataframe no need to run this
model=glm(Class.Label~.,family=binomial,data=train)#AIC: 10 
#here you will get this, no need to worry it means GLM model will fitted
#Warning message:
# glm.fit: fitted probabilities numerically 0 or 1 occurred 
summary(model)
model_pred=zapsmall( predict(model, newdata = train, type = "response") )
View(model_pred)
model_pred
test_pred = predict(model, newdata = test, type = "response")
test_pred=round(test_pred,1); test_pred
View(test_pred)
length(test_pred)
write.csv(test_pred,"Test.csv")
