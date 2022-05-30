library("coin")
library("ctree")
library(party)
setwd("/Users/klng 1/R/ML")
iris=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")

iris<- na.omit(iris) 
set.seed(1234)
#The random seed is set to a fixed value below to make the results reproducible. 
index <- sample(2, nrow(iris), replace=T, prob=c(0.7, 0.3)) # training sets [70%], test set [30%], CANNOT set replace=F
index
trainData <- iris[index==1,] 
head(trainData)
testData <- iris[index==2,] 
head(testData)
library(party)
myFormula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width 
iris_ctree <- ctree(myFormula, data=trainData) 
# check the prediction 
#table(predict(iris_ctree), trainData$Species) 
iris.pred<-table(predict(iris_ctree), trainData$Species)
iris.pred
acc<-round(sum(diag(iris.pred)) / sum(iris.pred)*100,2)
acc # signif( (40+37+31)/ (40+37+31+3+1),4)
tr <- treeresponse(iris_ctree, newdata = iris[1:3,])
tr

plot(iris_ctree)
plot(iris_ctree, type="simple")
# prediction using test data 
testPred <- predict(iris_ctree, newdata = testData) 
table(testPred, testData$Species) 
######## end of ctree ######################

######## using rparty
library(rpart)
data(iris)
index=ceiling(0.1*nrow(iris))
index
test.index=sample(1:nrow(iris),index)
iris.trainset=iris[-test.index,]  #use ‘-’ sign
iris.testset=iris[test.index,]
head(iris.trainset)
head(iris.testset)
#myFormula <- Species ~ .
iris_rpart <- rpart(Species ~ ., data=iris.trainset) ##later, we will do “Parameter optimization”
#iris_rpart <- rpart(Species ~ ., data=iris.trainset) , control=rpart.control(minsplit=5,cp=0.0001, maxdepth=30))
plot(iris_rpart)
text(iris_rpart)# add lebels
attributes(iris_rpart)
print(iris_rpart$frame) 

print(iris_rpart$cptable)
summary(iris_rpart)# three values of CP 
species.trainset=iris$Species[-test.index]
train.predict= factor(predict(iris_rpart,iris.trainset,type='class'),levels=levels(species.trainset))
table.trainset=table(species.trainset,train.predict)
table.trainset
acc.rpart<-round(sum(diag(table.trainset)) / sum(table.trainset)*100,2)
acc.rpart
species.testset=iris$Species[test.index]
test.predict= factor(predict(iris_rpart,iris.testset,type='class'),levels=levels(species.testset))
table.testset=table(species.testset,test.predict)
table.testset
acc.rpart.test<-round(sum(diag(table.testset)) / sum(table.testset)*100,2)
acc.rpart.test

############ rparty Parameter optimization 
?rpart.control
#minsplit = the minimum number of observations that must exist in a node in order for a split to be attempted
iris_rpart2<- rpart(Species ~ ., data=iris.trainset,control=rpart.control(minsplit=5,cp=0.0001, maxdepth=30))
plot(iris_rpart2)
text(iris_rpart2)
species.trainset=iris$Species[-test.index]
train.predict= factor(predict(iris_rpart,iris.trainset,type='class'),levels=levels(species.trainset))
table.trainset2=table(species.trainset,train.predict)
acc.rpart2<-round(sum(diag(table.trainset2)) / sum(table.trainset2)*100,2)
acc.rpart2
species.testset=iris$Species[test.index]
test.predict2= factor(predict(iris_rpart2,iris.testset,type='class'),levels=levels(species.testset))
table.testset2=table(species.testset,test.predict2)
table.testset2
acc.rpart.test2<-round(sum(diag(table.testset2)) / sum(table.testset2)*100,2)
acc.rpart.test2

######## end of rparty


