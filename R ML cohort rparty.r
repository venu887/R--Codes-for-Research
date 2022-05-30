# using DT to classify cancer cohort
setwd("/Users/klng 1/R/ML")
ov=read.csv("./OV.csv",header=TRUE) #CANNOT use header=TRUE, otherwise
dim(ov); col=dim(ov)[2]; col
col_minus1=col-1; col_minus1
ov[,7]; names(ov);
ov$CL
ov$CL<-as.factor(ov$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(ov$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
#rows 2:23 [22 samples] stages1+2, rows 24-301 stages3+4
# data imbalance problem # size of '0' >> size of '1'
dim(ov)[1]
# for different cohort change the num. 22 to something else
index_s12=ceiling(0.3*22) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.3*(dim(ov)[1]-22)) # data imbalance ==> need to sample individually
index_s12; index_s34

######### sampling, build the training set and test set #######################
set.seed(2307) 
test.index_s12=sample(1:22, index_s12); 
# for different cohort change the num. 22 to something else+1 
test.index_s34=sample(23:dim(ov)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)

trainset=ov[-test.index_s12,]; dim(trainset)  #use ‘-’ sign to remove entries
trainset=trainset[-test.index_s34,]; dim(trainset)
trainset = trainset[,-1]

testset12=ov$CL[test.index_s12]; length(testset12)
testset34=ov$CL[test.index_s34]; length(testset34)
testset=unlist(list(ov.testset12, testset34))
length(testset)

RPART.trainset=ov[-test.index_s12,]; dim(RPART.trainset)  #use ‘-’ sign to remove entries
RPART.trainset=RPART.trainset[-test.index_s34,]; dim(RPART.trainset)
RPART.trainset = RPART.trainset[,-1]

RPART.testset_s12=ov[test.index_s12,]; dim(RPART.testset_s12)
RPART.testset_s34=ov[test.index_s34,];  dim(RPART.testset_s34)
RPART.testset=rbind(RPART.testset_s12, RPART.testset_s34) #rbind to put together test data
dim(RPART.testset)
RPART.testset = RPART.testset[,-1]; dim(RPART.testset)

###  using DT 'party' for classification OR jump to rparty below ##############
######## using rparty
library(rpart)
myFormula <- CL~HOXA3.3200+HOXA5.3202+MEOX2.4223+PDGFD.80310+TCEAL7.56849 
#ov_rpart <- rpart(myFormula, data=DT.trainset) 
ov_rpart<- rpart(myFormula, data=trainset,control=rpart.control(minsplit=7,cp=0.00001, maxdepth=30))
##later, we will do “Parameter optimization”
#ov_rpart <- rpart(CL ~ ., data=ov.trainset) , control=rpart.control(minsplit=5,cp=0.0001, maxdepth=30))
plot(ov_rpart)
text(ov_rpart)# add lebels
attributes(ov_rpart)
print(ov_rpart$frame) 
summary(ov_rpart)# three values of CP 
#compute the rparty performance using training set
train.predict= factor(predict(ov_rpart, data=RPART.trainset,type='class'))
table.trainset=table(train.predict, RPART.trainset$CL)
table.trainset
acc.rpart.train<-round(sum(diag(table.trainset)) / sum(table.trainset)*100,2)
acc.rpart.train


#compute the rparty performance using test set
test.predict= factor(predict(ov_rpart,RPART.testset,type='class'))
#when use 'table' actual class label in the 1st arguemnt, prediction in the 2nd arguemnt, then FP, FN is correct
table.testset=table(testset, test.predict)
table.testset
acc.rpart.test<-signif(sum(diag(table.testset)) / sum(table.testset)*100,3)
acc.rpart.test

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
#out2 "0.923" "0.976" "0.286" "0.943" "0.959" "0.5" "0.164" "0.714" "0.0238"
