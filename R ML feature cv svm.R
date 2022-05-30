#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
setwd("/Users/klng 1/R/ML")
library(caret); set.seed(100)
iris=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")
#Pick the first two species for binary classification #Iris-setosa:50; Iris-versicolor:50  NOT USING Iris-virginica:50 
iris012=as.numeric(iris$Species)-1  #relabel the species = 0, 1 and 2, iris012
# need to subtract 1, else the label runs from 1 to 3 instead of 0 to 2.
IRIS=cbind(iris,iris012); dim(IRIS)
IRIS=IRIS[,-5]; dim(IRIS) # rm the $Species column 5 before using 'train '[character not num]
IRIS100=IRIS[1:100,]; dim(IRIS100)
IRIS100$iris012<-as.factor(IRIS100$iris012)
prop.table(table(IRIS100$iris012))
########################### cross-validation using svm ##########################
Control <- trainControl(method = "cv",number = 10)
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated
# SVMRadial
svm.m <- train(iris012 ~., data = IRIS100, method = "svmRadial", trControl = Control,  preProcess = c("center","scale"))
varImp(svm.m)
########### classification using SVM
setwd("/Users/klng 1/R/ML")
iris=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")
library(e1071)
iris100=iris[1:100,]; dim(iris100)
set.seed(23072307)
svm.index=ceiling(0.2*nrow(iris100))
svm.index
svm.test.index=sample(1:nrow(iris100), svm.index)
svm.trainset=iris100[-svm.test.index,]  #use ‘-’ sign
svm.testset=iris100[svm.test.index,]
#head(svm.trainset)   #head(svm.testset)
svm.model<-svm(Species ~ Petal.Width+Petal.Length, data=svm.trainset, type='C-classification',cost=10, gamma=10)
# gamma is a parameter for non linear hyperplanes. The higher the gamma value it tries to exactly 
#fit the training data set # C is the penalty parameter of the error term. It controls the trade off 
#between smooth decision boundary and classifying the training points correctly.
svm.pred<-predict(svm.model, svm.testset[,-5])
#svm.pred
table.svm.test=table(pred=svm.pred,true=svm.testset[,5])
table.svm.test
acc.svm<-round(sum(diag(table.svm.test)) / sum(table.svm.test)*100,2)
acc.svm # if you run the code again, one will get a different acc.svm value, which may be higher or lower
