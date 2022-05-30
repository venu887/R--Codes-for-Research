#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
setwd("J:/V/R")
rm(list = ls())
library(caret); set.seed(100)
ov=read.csv("ov miR_ML.csv",header=TRUE) #CANNOT use header=TRUE, otherwise
dim(ov); col=dim(ov)[2]; col
col_minus1=col-1; col_minus1
ov[,9]; names(ov);
ov$CL
#ov$CL<-as.factor(ov$CL)
fraction=signif(prop.table(table(ov$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
########################### cv using LR ##########################
Control <- trainControl(method = "cv",number = 10)
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated
# LR
# 1st column is TCGA ID, last columnis class label, skip these for cv ==> 2:col_minus1
logit.M <- train(ov[,2:col_minus1], ov[,col], method = "glmnet",
                 trControl = Control, tuneLength=2)
varImp(logit.M)
#####################################################
#rows 2:23 [22 samples] stages1+2, rows 24-301 stages3+4
# data imbalance problem # size of '0' >> size of '1'
dim(ov)[1]

index_s12=ceiling(0.3*72) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.3*(dim(ov)[1]-72)) # data imbalance ==> need to sample individually
index_s12; index_s34
######### sampling, build the training set and test set #######################
set.seed(2307) 
test.index_s12=sample(1:72, index_s12); 
test.index_s34=sample(73:dim(ov)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)

LR.trainset=ov[-test.index_s12,]; dim(LR.trainset)  #use ‘-’ sign to remove entries
LR.trainset=LR.trainset[-test.index_s34,]; dim(LR.trainset)
LR.trainset = LR.trainset[,-1]

LR.testset_s12=ov[test.index_s12,]; dim(LR.testset_s12)
LR.testset_s34=ov[test.index_s34,];  dim(LR.testset_s34)
LR.testset=rbind(LR.testset_s12, LR.testset_s34) #rbind to put together test data
dim(LR.testset)
#LR.testset = LR.testset[,-1]; dim(LR.testset)
names(ov)
myFormula_ov <- CL~hsa.mir.1911+hsa.mir.449a+hsa.mir.605+hsa.mir.449b+hsa.mir.34b+hsa.mir.34c+hsa.mir.767
#myFormula_ucec<-CL~hsa.mir.891b+hsa.mir.137+hsa.mir.1197+hsa.mir.675+hsa.mir.211+hsa.mir.548y+hsa.mir.1251+hsa.mir.122+hsa.mir.490+hsa.mir.599+hsa.mir.184+hsa.mir.889+hsa.mir.483+hsa.mir.323b
#myFormula_cesc<-CL~hsa.mir.144+hsa.mir.486+hsa.mir.508+hsa.mir.526b+hsa.mir.522+hsa.mir.935+hsa.mir.204+hsa.mir.192+hsa.mir.514b+hsa.mir.506

model=glm(CL~hsa.mir.1911+hsa.mir.449a+hsa.mir.605+hsa.mir.449b+hsa.mir.34b+hsa.mir.34c+hsa.mir.767, family=binomial,data=LR.trainset)
summary(model) 

names(LR.trainset) #when use LR do not consider 1st column
model_pred=zapsmall( predict(model, newdata =LR.trainset , type = "response") ); model_pred
test_pred = predict(model, newdata = LR.testset, type = "response")
test_pred=round(test_pred,1); test_pred

length(test_pred)
head(test_pred)
table.LR.test=table(pred=test_pred,true=LR.testset$CL)
table.LR.test
str(table.LR.test)
acc.LR<-round(sum(diag(table.LR.test)) / sum(table.LR.test)*100,2)
acc.LR
TN=table.LR.test[1,1]; TN
TP=table.LR.test[2,2]; TP
FN=table.LR.test[2,1]; FN
FP=table.LR.test[1,2]; FP
#accuracy, Q, sensitivity, SN, specificity, SP, F1 measure, F1, positive predictive rate, 
#STP , negative predictive rate, STN , false positive rate (alpha), and 
#false negative rate (beta), they are defined as 
Q = (TP+TN)/(TP+TN+FP+FN); Q
SN = TP/(TP+FN); SN
SP = TN/(TN+FP); SP
STP = TP/(TP+FP); STP
F1 =1/( 1/2*(1/SN + 1/STP) ); F1
#negative predictive rate,
STN = TN/(TN+FN); STN
MCC = ((TP+TN)-(FP*FN))/ ( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ); MCC
#False positive rate (alpha),
alpha= FP / (FP + TN) #= 1 – specificity
alpha
#False negative rate (beta),
beta= FN / (TP + FN); beta 
