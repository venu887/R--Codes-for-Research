rm(list = ls())
setwd("J:/V/R")
iris=read.csv("cesc miR_logcounts-ML.csv",header=TRUE) #CANNOT use header=TRUE
#names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")
dim(iris)
library(e1071)
#iris100=iris[1:100,]; dim(iris100)
set.seed(23072307)
svm.index=ceiling(0.2*nrow(iris))
svm.index
svm.test.index=sample(1:nrow(iris), svm.index)
sort(svm.test.index)
svm.trainset=iris[-svm.test.index,]  #use ‘-’ sign
svm.testset=iris[svm.test.index,]
dim(svm.testset)
dim(svm.trainset)
names(iris)
#myFormula_ov <- CL~hsa.mir.1911+hsa.mir.449a+hsa.mir.605+hsa.mir.449b+hsa.mir.34b+hsa.mir.34c+hsa.mir.767
#myFormula_ucec<-CL~hsa.mir.891b+hsa.mir.137+hsa.mir.1197+hsa.mir.675+hsa.mir.211+hsa.mir.548y+hsa.mir.1251+hsa.mir.122+hsa.mir.490+hsa.mir.599+hsa.mir.184+hsa.mir.889+hsa.mir.483+hsa.mir.323b
myFormula_cesc<-CL~hsa.mir.144+hsa.mir.486+hsa.mir.508+hsa.mir.526b+hsa.mir.522+hsa.mir.935+hsa.mir.204+hsa.mir.192+hsa.mir.514.1+hsa.mir.509.2+hsa.mir.509.3
svm.model<-svm(myFormula_cesc, data=svm.trainset, type='C-classification',cost=10, gamma=10)
# gamma is a parameter for non linear hyperplanes. The higher the gamma value it tries to exactly 
#fit the training data set # C is the penalty parameter of the error term. It controls the trade off 
#between smooth decision boundary and classifying the training points correctly.
svm.pred<-predict(svm.model, svm.testset[,-13]) #change from 5 to 7
#svm.pred
table.svm.test=table(true=svm.testset$CL, pred=svm.pred) #change from 5 to 7
table.svm.test
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm # if you run the code again, one will get a different acc.svm value, which may be higher or lower

tuned<-tune.svm(CL ~ ., data=svm.trainset[,2:13], gamma=10^(-8:1), cost=2^(0:4)) # differ by 3 orders of magnitude
summary(tuned)
plot(tuned, data = ov)
svm.model2<-svm(CL ~ ., data=svm.trainset[,2:13], kernal="radial", gamma=10,cost=10)
svm.pred2<-predict(svm.model2, svm.testset[,-1])

table.svm.test2=table( pred=svm.pred2, true=svm.testset$CL)
table.svm.test2
str(table.svm.test2)
#miss classification rate 
1- sum(diag(table.svm.test2))/sum(table.svm.test2)
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm # if you run the code again, one will get a different acc.svm value, which may be higher or lower
TN=table.svm.test2[1,1]; TN
TP=table.svm.test2[2,2]; TP
FN=table.svm.test2[2,1]; FN
FP=table.svm.test2[1,2]; FP
#accuracy, Q, sensitivity, SN, specificity, SP, F1 measure, F1, positive predictive rate, 
#STP , negative predictive rate, STN , false positive rate (alpha), and 
#false negative rate (beta), they are defined as 
Q = (TP+TN)/(TP+TN+FP+FN); signif(Q,3)
SN = TP/(TP+FN); signif(SN,3)
SP = TN/(TN+FP); signif(SP,3)
STP = TP/(TP+FP); signif(STP,3)
F1 =1/( 1/2*(1/SN + 1/STP) ); F1
#negative predictive rate,
STN = TN/(TN+FN); STN
MCC = ((TP+TN)-(FP*FN))/ ( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ); MCC
#False positive rate (alpha),
alpha= FP / (FP + TN) #= 1 - specificity
alpha
#False negative rate (beta),
beta= FN / (TP + FN); beta 
set.seed(123)

