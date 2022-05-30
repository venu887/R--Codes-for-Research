# using adaboost to classify cancer cohort
#https://www.datatechnotes.com/2018/03/classification-with-adaboost-model-in-r.html 
# rfSRC, https://cran.r-project.org/web/packages/randomForestSRC/index.html
# use LASSO 
#https://stats.stackexchange.com/questions/72251/an-example-lasso-regression-using-glmnet-for-binary-outcome
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r 
setwd("/Users/klng 1/R/ML")
library(lattice)
library(gtable)
library(munsell)
library(scales)
library(colorspace)
library(lazyeval)
library(ggplot2)
library(caret)
library(adabag)  
set.seed(100)

abc=read.csv("./abc_demo2.csv",header=TRUE,row.names = 1) #CANNOT use header=TRUE, otherwise
dim(abc); col=dim(abc)[2]; col
#col_minus1=col-1; col_minus1
abc[,6]; names(abc);
abc$CL
abc$CL<-as.factor(abc$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(abc$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'

#####################################################
# data imbalance problem # size of '0' >> size of '1'
dim(abc)[1]
index_s12=ceiling(0.4*22) #use 0.2 ==> remove test set index from raw dataset
index_s34=ceiling(0.4*(dim(abc)[1]-22)) # data imbalance ==> need to sample individually
index_s12; index_s34
######### sampling, build the training set and test set #######################
set.seed(2307) 
test.index_s12=sample(1:22, index_s12); 
test.index_s34=sample(23:dim(abc)[1], index_s34); 
sort(test.index_s12); sort(test.index_s34); length(test.index_s12); length(test.index_s34)

RF.trainset=abc[-test.index_s12,]; dim(RF.trainset)  #use ‘-’ sign to remove entries
RF.trainset=RF.trainset[-test.index_s34,]; dim(RF.trainset)


RF.testset_s12=abc[test.index_s12,]; dim(RF.testset_s12)
RF.testset_s34=abc[test.index_s34,];  dim(RF.testset_s34)
RF.testset=rbind(RF.testset_s12, RF.testset_s34) #rbind to put together test data
dim(RF.testset)

# use adaboost
myFormula <- CL~ HOXA3.3200+HOXA5.3202+MEOX2.4223+PDGFD.80310+TCEAL7.56849 
train.adaboost<-boosting(myFormula, data=RF.trainset, boos=TRUE, mfinal=5)
# iterate mfinal times. with mfinal > 5 does not imply better acc BUT test it
ada.pred=predict.boosting(train.adaboost,newdata=RF.testset)
ada.pred$confusion
ada.pred$error
ada.acc=1-ada.pred$error
ada.acc

confusion.matrix=table(actual=RF.testset$CL, pred= ada.pred$class)
confusion.matrix

TN=confusion.matrix[1,1]; TN
TP=confusion.matrix[2,2]; TP
FN=confusion.matrix[2,1]; FN
FP=confusion.matrix[1,2]; FP
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

