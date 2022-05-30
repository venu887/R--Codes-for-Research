setwd("/Users/user/Documents/R/ML")
iris=read.csv("./iris.data", header=F)
dim(iris)
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")
head(iris,3);
tail(is.na(iris))
iris012=as.numeric(as.factor(iris$Species))-1
iris012
IRIS=cbind(iris,iris012) # insert the ‘class label’ column
dim(IRIS)

iris100=IRIS[1:100, -5] # delete the 'Species' column [charaters]
head(iris100); dim(iris100)
#cov=cov(iris100); cov
#eigen(cov)
#prcomp(cov) # no output eigenvalues but eignevectors

set.seed(2323) 
indexes<-sample(nrow(iris100),0.8*nrow(iris100),replace = F) # training set 80%  100*0.8=80 entries, test set 20%  20 entries 
indexes
train<-iris100[indexes,]; dim(train)
test<-iris100[-indexes,]; dim(test)

myFormula <- iris012 ~ .
svm.model<-svm(myFormula, data=train, type='C-classification',cost=10, gamma=10)
svm.pred<-predict(svm.model, test[,-5]) #no need to use iris012 column
#svm.pred
table.svm.test=table(pred=svm.pred,true=test[,5]) #change from 5 to 7
table.svm.test
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm #


