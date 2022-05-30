#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
rm(list = ls())
setwd("J:/V/R")

OV_20<-read.csv("ov miR_ML.csv", row.names = 1)# CSV file Raw data or statewide information of RSEM values
acc<-read.csv("OV 20 DEG.csv", row.names = 1) # CSV file after limma results DEGs
UPdeg<-acc[acc$logFC>= 1 & acc$adj.P.Val<0.05,] #UPREGULATED GENES
Downdeg<-acc[acc$logFC<= -1 & acc$adj.P.Val<0.05,] #DOWNREGULATED GENES
acc_final<-rbind(UPdeg,Downdeg)
rownames(acc_final)
library(stringr)

b<-merge(OV_20,acc_final, by = "row.names", all.x = FALSE,all.y = FALSE)
b[,1] # after merrging we can see row names is in first column 
rownames(b)<-b[,1] # we have to replace row names 
b[,1]
View(b)
b<-b[,-1]
new_b<-b[,1:300]
bt<-t(new_b)
write.csv(rownames(bt), "Class Label.csv")
# arrange in excel change patient names as row names 
# add new column Class Label S12 as negative "0" S34 as Positive "1"
CL<-read.csv("Class Label.csv", row.names = 1)
View(CL)
class(CL)
d<-cbind(bt,CL)
colnames(d)
write.csv(d, "OV for svm.csv")

#install.packages("caret")
library(caret); set.seed(100)
ov=read.csv("Cancer RSEM.csv",header=TRUE) #CANNOT use header=TRUE, otherwise
#ggplot 
library(ggplot2)
names(ov)
qplot(hsa.mir.1911,data = ov, color=CL)


dim(ov); col=dim(ov)[2]
ov[,13]; names(ov)
ov$CL
any(is.na(ov))
ov$CL<-as.factor(ov$CL)
signif(prop.table(table(ov$CL)),3) #rows 1-168 stages1+2, rows 169-176 stages3+4
# data imbalance problem # size of '0' >> size of '1'
# 1st column is TCGA ID, skip that for cv
########################### cross-validation using svm ##########################
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated # SVMRadial
Control <- trainControl(method = "cv",number = 10) # cross validation
svm.m <- train(CL ~., data = ov[,2:dim(ov)[2]], method = "svmRadial", trControl = Control,  preProcess = c("center","scale"))
varImp(svm.m)
plot(varImp(svm.m))
summary(svm.m)
########### classification using SVM
#install.packages("e1071")
library(e1071)
set.seed(123)
index_S12=ceiling(0.3*233) #use 0.2 ==> remove test set index from raw dataset
index_S34=ceiling(0.3*(dim(ov)[1]-233)) # data imbalance ==> need to sample individually
index_S12; index_S34
#extract selected index data from s12 and s34
test.index_S12=sample(1:233, index_S12)
test.index_S34=sample(234:dim(ov)[1], index_S34); 
sort(test.index_S12); sort(test.index_S34); length(test.index_S12); length(test.index_S34)
# isolate test set data 
svm.testset_S12=ov[test.index_S12,]; dim(svm.testset_S12)
svm.testset_S34=ov[test.index_S34,];  dim(svm.testset_S34)
svm.testset=rbind(svm.testset_S12, svm.testset_S34) #rbind to put together test data
# isolate train set data by removing test set indexes
svm.trainset1=ov[-test.index_S12,]; dim(svm.trainset1) #-ve class  #use '-' sign to remove entries
svm.trainset=svm.trainset1[-test.index_S34,]; dim(svm.trainset) #+ve class 
summary(svm.trainset)
names(svm.trainset) #when use SVM do not consider 1st column
dim(svm.trainset)
dim(svm.testset1)
#head(svm.trainset)   #head(svm.testset)
svm.model<-svm(CL~. , data=svm.trainset[,2:9], type='C-classification')
summary(svm.model)
svmfit<-svm(CL~.,data=svm.trainset[,2:9], kernel="linear", scale=FALSE)
summary(svmfit)
print(svmfit)
plot(svmfit, data = svm.trainset[,2:9], hsa.mir.605~hsa.mir.767)

# gamma is a parameter for non linear hyperplanes. The higher the gamma value it tries to exactly 
#fit the training data set # C is the penalty parameter of the error term. It controls the trade off 
#between smooth decision boundary and classifying the training points correctly.
svm.pred<-predict(svm.model, svm.testset1[,-1])
#svm.pred
table.svm.test=table(pred=svm.pred,true=svm.testset1$CL)
table.svm.test
#miss classification rate 
1- sum(diag(table.svm.test))/sum(table.svm.test)
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm # if you run the code again, one will get a different acc.svm value, which may be higher or lower
########## tuning svm parameters ################
#hyper parameter optimization to predict the best kernel model 
set.seed(123)
#t<-tune(svm, CL~., data = ov, ranges = list(epsilon= seq(0, 1, 0.1)), cost=10^(2:10))
#plot(t)
tuned<-tune.svm(CL ~ ., data=svm.trainset[,2:9], gamma=10^(-3:1), cost=10^(-2:2)) # differ by 3 orders of magnitude
summary(tuned)
plot(tuned, data = ov)
svm.model2<-svm(CL ~ ., data=svm.trainset[,2:9], kernal="radial", gamma=0.01,cost=10)
svm.pred2<-predict(svm.model2, svm.testset1[,-1])
head(svm.pred2)
table.svm.test2=table(pred=svm.pred2,true=svm.testset1$CL)
table.svm.test2
str(table.svm.test2)
#miss classification rate 
1- sum(diag(table.svm.test2))/sum(table.svm.test2)
acc.svm2<-round(sum(diag(table.svm.test2)) / sum(table.svm.test2)*100,2)
acc.svm2
TN=table.svm.test2[1,1]; TN
TP=table.svm.test2[2,2]; TP
FN=table.svm.test2[2,1]; FN
FP=table.svm.test2[1,2]; FP
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
alpha= FP / (FP + TN) #= 1 - specificity
alpha
#False negative rate (beta),
beta= FN / (TP + FN); beta 