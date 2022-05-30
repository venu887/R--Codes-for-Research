# using DT to classify cancer cohort
rm(list = ls())
setwd("J:/V/R")
library(caret); set.seed(100)
xyz=read.csv("cesc miR_ML.csv",header=TRUE) #CANNOT use header=TRUE, otherwise
dim(xyz); col=dim(xyz)[2]; col
col_minus1=col-1; col_minus1
xyz[,12]; names(xyz);
xyz$CL
xyz$CL<-as.factor(xyz$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(xyz$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
########################### no cv using DT ##########################
#Control <- trainControl(method = "cv",number = 10)
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated
# DT
# 1st column is TCGA ID, last columnis class label, skip these for cv ==> 2:col_minus1
#logit.M <- train(xyz[,2:col_minus1], xyz[,col], method = "glmnet",
#                 trControl = Control, tuneLength=2)
#varImp(logit.M)
#####################################################
#rows 2:23 [22 samples] stages1+2, rows 24-301 stages3+4
# data imbalance problem # size of '0' >> size of '1'
dim(xyz)[1]
# for different cohort change the num. 22 to something else
index_s12=ceiling(0.2*233) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.2*(dim(xyz)[1]-233)) # data imbalance ==> need to sample individually
index_s12; index_s34

######### sampling, build the training set and test set #######################
set.seed(2307) 
test.index_s12=sample(1:233, index_s12); 
# for different cohort change the num. 22 to something else+1 
test.index_s34=sample(234:dim(xyz)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)

DT.trainset=xyz[-test.index_s12,]; dim(DT.trainset)  #use ‘-’ sign to remove entries
DT.trainset=DT.trainset[-test.index_s34,]; dim(DT.trainset)
DT.trainset = DT.trainset[,-1]

DT.testset_s12=xyz[test.index_s12,]; dim(DT.testset_s12)
DT.testset_s34=xyz[test.index_s34,];  dim(DT.testset_s34)
DT.testset=rbind(DT.testset_s12, DT.testset_s34) #rbind to put together test data
dim(DT.testset)
DT.testset = DT.testset[,-1]; dim(DT.testset)
###  using DT 'party' for classification  ##############
library("coin")
library(party)
#myFormula_ov <- CL~hsa.mir.1911+hsa.mir.449a+hsa.mir.605+hsa.mir.449b+hsa.mir.34b+hsa.mir.34c+hsa.mir.767
#myFormula_ucec<-CL~hsa.mir.891b+hsa.mir.137+hsa.mir.1197+hsa.mir.675+hsa.mir.211+hsa.mir.548y+hsa.mir.1251+hsa.mir.122+hsa.mir.490+hsa.mir.599+hsa.mir.184+hsa.mir.889+hsa.mir.483+hsa.mir.323b
#myFormula_cesc<-CL~hsa.mir.144+hsa.mir.486+hsa.mir.508+hsa.mir.526b+hsa.mir.522+hsa.mir.935+hsa.mir.204+hsa.mir.192+hsa.mir.514b+hsa.mir.506
xyz.ctree <- ctree(myFormula_cesc, data=DT.trainset) # xyz.ctree = MODEL
# check the prediction using training set
table.trainset<-table(true=DT.trainset$CL, pred=predict(xyz.ctree))
table.trainset
plot(xyz.ctree)
plot(xyz.ctree, type="simple")
######## end of ctree ######################
#compute the DT performance using training set
#train.predict= factor(predict(xyz.ctree,DT.trainset,type='response'))
#when use 'table' actual class label in the 1st arguemnt, prediction in the 2nd arguemnt, then FP, FN is correct
#table.trainset=table(true=DT.trainset$CL, pred=train.predict)
#table.trainset

TN=table.trainset[1,1]; TN
TP=table.trainset[2,2]; TP
FN=table.trainset[2,1]; FN
FP=table.trainset[1,2]; FP
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

##### compute the DT performance using test set
#test_pred <- predict(xyz.ctree, newdata = DT.testset) 
#table.testset<-table(true=DT.testset$CL, pred=predict(xyz.ctree))
test.predict= factor(predict(xyz.ctree,DT.testset,type='response'))
#when use 'table' actual class label in the 1st arguemnt, prediction in the 2nd arguemnt, then FP, FN is correct
table.testset=table(true=DT.testset$CL, pred=test.predict)
table.testset # there is NO '0' class prediction, table has '1' prediction !!
acc.test<-signif(sum(diag(table.testset)) / sum(table.testset)*100,3)
acc.test
table(test_pred, DT.testset$CL)
TN=table.testset[1,1]; TN
TP=table.testset[2,2]; TP
FN=table.testset[2,1]; FN
FP=table.testset[1,2]; FP
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

#[,1]    [,2]    [,3] [,4]    [,5]    [,6]  [,7]   [,8]    [,9]    
#out1 "Q"     "SN"    "SP"    "STP"   "F1"    "STN" "MCC"  "alpha" "beta"  
#out2 "0.901" "0.976" "0"     "0.921" "0.948" "0"   "0.21" "1"      "0.0238"
#after optimized the parameters, below is what I got
#minsplit=5,cp=0.0001, maxdepth=30
#out2 "0.923" "0.976" "0.286" "0.943" "0.959" "0.5" "0.164" "0.714" "0.0238" 
#minsplit=5,cp=0.00001, maxdepth=10
#out2 "0.923" "0.976" "0.286" "0.943" "0.959" "0.5" "0.164" "0.714" "0.0238"
#minsplit=5,cp=0.00001, maxdepth=30
#out2 "0.923" "0.976" "0.286" "0.943" "0.959" "0.5" "0.164" "0.714" "0.0238"
#minsplit=7,cp=0.00001, maxdepth=30
#out2 "0.923" "0.976" "0.286" "0.943" "0.959" "0.5" "0.164" "0.714" "0.0238"
