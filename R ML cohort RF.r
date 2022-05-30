# using DT to classify cancer cohort
rm(list = ls())
setwd("J:/V/R")
library(randomForest)
set.seed(100)
v=read.csv("Cancer RSEM.csv") #CANNOT use header=TRUE, otherwise
dim(v); col=dim(v)[2]; col
col_minus1=col-1; col_minus1
v[,17]; names(v);
v$CL
v$CL<-as.factor(v$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(v$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
########################### no cv using DT ##########################
#Control <- trainControl(method = "cv",number = 10)
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated
# DT
# 1st column is TCGA ID, last columnis class label, skip these for cv ==> 2:col_minus1
#logit.M <- train(v[,2:col_minus1], v[,col], method = "glmnet",trControl = Control, tuneLength=2)
#varImp(logit.M)
#####################################################
# data imbalance problem # size of '0' >> size of '1'
#dim(v)[1]
#index_s12=ceiling(0.2*233) #use 0.2 ==> remve test set index from raw dataset
#index_s34=ceiling(0.2*(dim(v)[1]-233)) # data imbalance ==> need to sample individually
#index_s12; index_s34
######### sampling, build the training set and test set #######################
set.seed(23072307)
rf.index=ceiling(0.4*nrow(v))
rf.index
rf.test.index=sample(1:nrow(v), rf.index)
rf.trainset=v[-rf.test.index,]  #use '-' sign
rf.testset=v[rf.test.index,]
names(v)
myFormula <- CL~ENSG00000140465.12+ENSG00000140505.6+ENSG00000153802.10+ENSG00000169594.11+ENSG00000170454.5+ENSG00000181433.8+ENSG00000184330.10+ENSG00000198883.10+ENSG00000213231.11+ENSG00000223850.1+ENSG00000227121.1+ENSG00000233639.4+ENSG00000250920.1+ENSG00000260676.4+ENSG00000268864.3
names(v)
rf <- randomForest(myFormula,data=rf.trainset,ntree=500,nodesize=3,proximity=TRUE,corr.bias=TRUE) 
print(rf) 
attributes(rf) 
plot(rf) 
importance(rf) 
varImpPlot(rf)
# apply rf model on test set#rf.pred <- predict(rf, newdata=RF.testset) 
#table(rf.pred, RF.testset$CL) #plot(margin(rf, RF.testset$CL)) 
# compute the performance of rf on training set 
rf.pred = predict(rf, data=rf.trainset)
table.trainset.rf<-table(rf.pred,rf.trainset$CL)
table.trainset.rf
acc.rf<-signif(sum(diag(table.trainset.rf)) / sum(table.trainset.rf)*100,3)
acc.rf
library(pROC)
rf_roc<-roc(RF.trainset$CL,rf$votes[,2],plot=TRUE,print.auc=TRUE)
plot(rf_roc, col=2, main="CESC RF", asp=NA)
#par(pty = "s")
abline(h=0.71, v=0.45)
legend("bottomright",0.2, 0.2,  legend=c("AUC=0.6509"), lwd=1)
legend(0., 0.2, c("AUC=0.6509"),1:3)
auc(rf.roc)

test.predict= factor(predict(rf,rf.testset,type='class'),levels=levels(rf.testset))
#when use 'table'  actual class label in the 1st arguemnt, prediction in the 2nd arguemnt, then FP, FN is correct
table.testset=table(rf.testset,test.predict)
table.testset
acc.rf.test<-round(sum(diag(table.testset)) / sum(table.testset)*100,2)
acc.rf.test
roc.test <- roc(RF.trainset$CL,rf$votes[,2])
plot(rf.test)
auc(roc.test)

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

#[,1]    [,2] [,3]    [,4]    [,5]    [,6]  [,7]    [,8]    [,9]  
#out1 "Q"     "SN" "SP"    "STP"   "F1"    "STN" "MCC"   "alpha" "beta"
#out2 "0.934" "1"  "0.143" "0.933" "0.966" "1"   "0.369" "0.857" "0"   

#after optimized the parameters, below is what I got
#minsplit=5,cp=0.0001, maxdepth=30