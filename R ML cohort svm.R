#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
setwd("J:/V/R")
xy=read.csv("cesc miR_logcounts-ML.csv",header=TRUE) #CANNOT use header=TRUE
library(caret)
library(ggplot2)
library(e1071); set.seed(100)
dim(xy); col=dim(xy)[2]; col
xy[,col]; names(xy)  #
xy$CL
xy$CL<-as.factor(xy$CL)
imbalance=signif(prop.table(table(xy$CL)),3) #rows 1-168 stages1+2, rows 169-176 stages3+4
# data imbalance problem # size of '0' >> size of '1'
r_neg_over_pos= imbalance[1]/imbalance[2]; signif(r_neg_over_pos*100,3)
# 1st column is TCGA ID, skip that for cv
########################### cross-validation using svm ##########################
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated # SVMRadial
Control <- trainControl(method = "cv",number = 10)
svm.m <- train(CL ~., data = xy[,2:dim(xy)[2]], method = "svmRadial", trControl = Control,  preProcess = c("center","scale"))
varImp(svm.m)
########### classification using SVM
# SVM package
set.seed(2323)
# for the present cohort there are 22 records for stage12
# you need to open the file to find out 22, R code does not do that
index_s12=ceiling(0.2*233) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.2*(dim(xy)[1]-233)) # data imbalance ==> need to sample individually
index_s12; index_s34
# need to change sample from 1:22, and start from 23:dim(xy)
test.index_s12=sample(1:233, index_s12); 
test.index_s34=sample(234:dim(xy)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)
svm.trainset=xy[-test.index_s12,]  #use ‘-’ sign to remove entries
svm.trainset=svm.trainset[-test.index_s34,]; dim(svm.trainset)

names(svm.trainset) #when use SVM do not consider 1st column
svm.testset_s12=xy[test.index_s12,]; dim(svm.testset_s12)
svm.testset_s34=xy[test.index_s34,];  dim(svm.testset_s34)
svm.testset=rbind(svm.testset_s12, svm.testset_s34) #rbind to put together test data
dim(svm.testset)
#head(svm.trainset)   #head(svm.testset)
#myFormula <- CL~TCEAL7|56849+PDGFD|80310+HOXA5|3202+MEOX2|4223+HOXA3|3200
#Change it to 
myFormula_cesc<-CL~hsa.mir.144+hsa.mir.486+hsa.mir.508+hsa.mir.526b+hsa.mir.522+hsa.mir.935+hsa.mir.204+hsa.mir.192+hsa.mir.514.1+hsa.mir.509.2+hsa.mir.509.3
svm.model<-svm(myFormula_cesc , data=svm.trainset[,2:col], type='C-classification',cost=10, gamma=10)
svm.model
svmfit<-svm(myFormula_cesc,data=svm.trainset[,2:col], kernel="linear", scale=FALSE)
xy[,2:dim(xy)[2]]
#plot(svm.model, svm.trainset[,2:col], hsa.mir.1911~hsa.mir.449a)
# gamma is a parameter for non linear hyperplanes. The higher the gamma value it tries to exactly 
#fit the training data set # C is the penalty parameter of the error term. It controls the trade off 
#between smooth decision boundary and classifying the training points correctly.
svm.pred<-predict(svm.model, svm.testset[,-1])
#svm.pred
table.svm.test=table(true=svm.testset$CL,pred=svm.pred)
table.svm.test
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm # if you run the code again, one will get a different acc.svm value, which may be higher or lower
TN=table.svm.test[1,1]; TN
TP=table.svm.test[2,2]; TP
FN=table.svm.test[2,1]; FN
FP=table.svm.test[1,2]; FP
#accuracy, Q, sensitivity, SN, specificity, SP, F1 measure, F1, positive predictive rate, 
#STP , negative predictive rate, STN , false positive rate (alpha), and 
#false negative rate (beta), they are defined as 
Q = (TP+TN)/(TP+TN+FP+FN); signif(Q,3)
SN = TP/(TP+FN); signif(SN,3)
SP = TN/(TN+FP); signif(SP,3)
STP = TP/(TP+FP); signif(STP,3)
F1 =1/( 1/2*(1/SN + 1/STP) ); signif(F1,3)
#negative predictive rate,
STN = TN/(TN+FN); signif(STN,3)
MCC = ((TP+TN)-(FP*FN))/ ( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ); signif(MCC,3)
#False positive rate (alpha),
alpha= FP / (FP + TN) #= 1 – specificity
signif(alpha,3)
#False negative rate (beta),
beta= FN / (TP + FN); signif(beta,3)

out1=c("Q","SN","SP", "STP", "F1", "STN","MCC","alpha","beta")
out2=c(signif(Q,3),signif(SN,3),signif(SP,3), signif(STP,3),signif(F1,3),signif(STN,3),signif(MCC,3),signif(alpha,3),signif(beta,3))
out=rbind(out1, out2); out
########## tuning svm parameters ################
tuned<-tune.svm(CL ~ ., data=svm.trainset[,2:13], gamma=10^(-3:1), cost=10^(-2:2)) # differ by 3 orders of magnitude
summary(tuned)
svm.model2<-svm(CL ~ ., data=svm.trainset[,2:13], kernal="radial", gamma=0.001,cost=0.01)
svm.pred2<-predict(svm.model2, svm.testset[,-1])
head(svm.pred2)
table.svm.test2=table(pred=svm.pred2,true=svm.testset$CL)
table.svm.test2
str(table.svm.test2)
acc.svm2<-round(sum(diag(table.svm.test2)) / sum(table.svm.test2)*100,2)
acc.svm2
Q = (TP+TN)/(TP+TN+FP+FN); signif(Q,3)
SN = TP/(TP+FN); signif(SN,3)
SP = TN/(TN+FP); signif(SP,3)
STP = TP/(TP+FP); signif(STP,3)
F1 =1/( 1/2*(1/SN + 1/STP) ); signif(F1,3)
#negative predictive rate,
STN = TN/(TN+FN); signif(STN,3)
MCC = ((TP+TN)-(FP*FN))/ ( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ); signif(MCC,3)
#False positive rate (alpha),
alpha= FP / (FP + TN) #= 1 – specificity
signif(alpha,3)
#False negative rate (beta),
beta= FN / (TP + FN); signif(beta,3)

out1=c("Q","SN","SP", "STP", "F1", "STN","MCC","alpha","beta")
out2=c(signif(Q,3),signif(SN,3),signif(SP,3), signif(STP,3),signif(F1,3),signif(STN,3),signif(MCC,3),signif(alpha,3),signif(beta,3))
out=rbind(out1, out2); out

#kernal="radial", gamma=0.10,cost=10
#out1 "Q"     "SN"    "SP"     "STP"   "F1"     "STN"   "MCC"    "alpha"  "beta"
#out2 "0.986" "0.8"   "1"       "1"    "0.889"  "0.985" "0.235"  "0"      "0.2" 
# kernal="radial", gamma=0.01,cost=100
#out2 "0.958" "0.667" "0.971"  "0.5"  "0.571"   "0.985" "0.282"  "0.029"  "0.333"


