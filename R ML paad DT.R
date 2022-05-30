setwd("/Users/klng 1/R/ML")
paad=read.csv("cesc miR_logcounts-ML.csv",header=TRUE) #CANNOT use header=TRUE, otherwise
dim(paad); col=dim(paad)[2]; col
paad[,13]; names(paad)
paad$CL
paad$CL<-as.factor(paad$CL)
signif(prop.table(table(paad$CL)),3) #rows 1-168 stages1+2, rows 169-176 stages3+4
# data imbalance problem # size of '0' >> size of '1'
# 1st column is TCGA ID, skip that for cv

########### classification using SVM
library("coin")
#library("ctree")
library(party)
################# building training and test sets   #######################
set.seed(2307)
index_s12=ceiling(0.2*233) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.2*(dim(paad)[1]-233)) # data imbalance ==> need to sample individually
index_s12; index_s34
test.index_s12=sample(1:233, index_s12); 
test.index_s34=sample(234:dim(paad)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)
dt.trainset=paad[-test.index_s12,]  #use ‘-’ sign to remove entries
dt.trainset=dt.trainset[-test.index_s34,]; dim(dt.trainset)
names(dt.trainset) #when use SVM do not consider 1st column
dt.testset_s12=paad[test.index_s12,]; dim(dt.testset_s12)
dt.testset_s34=paad[test.index_s34,];  dim(dt.testset_s34)
dt.testset=rbind(dt.testset_s12, dt.testset_s34) #rbind to put together test data
dim(dt.testset)
################# building training and test sets   #######################
myFormula <- CL ~ . 
dt_ctree <- ctree(myFormula, data=dt.trainset[,2:13]) 
# check the prediction 
dt.pred<-table(predict(dt_ctree), dt.trainset$CL)
dt.pred
acc.train<-signif(sum(diag(dt.pred)) / sum(dt.pred)*100,3); acc.train
#tr <- treeresponse(dt_ctree, newdata = paad[1:3,]); tr
plot(dt_ctree)
plot(dt_ctree, type="simple")
dt.pred.test <- predict(dt_ctree, newdata = dt.testset) 
table.dt.test=table(pred=dt.pred.test,true=dt.testset$CL)
table.dt.test
acc.test<-signif(sum(diag(table.dt.test)) / sum(table.dt.test)*100,3)
acc.test # if you run the code again, one will get a different acc.svm value, which may be higher or lower
##########  DT performance ################
TN=table.dt.test[1,1]; TN
TP=table.dt.test[2,2]; TP
FN=table.dt.test[2,1]; FN
FP=table.dt.test[1,2]; FP
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

