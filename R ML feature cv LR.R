#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
setwd("J:/V/R")
library(caret); set.seed(100)
iris #CANNOT use header=TRUE, otherwise
#Pick the first two species for binary classification #Iris-setosa:50; Iris-versicolor:50  NOT USING Iris-virginica:50 
iris012=as.numeric(iris$Species)-1  #relabel the species = 0, 1 and 2, iris012
# need to subtract 1, else the label runs from 1 to 3 instead of 0 to 2.
IRIS=cbind(iris,iris012); dim(IRIS)
IRIS=IRIS[,-5]; dim(IRIS) # rm the $Species column 5 before using 'train '[character not num]
IRIS100=IRIS[1:100,]; dim(IRIS100)
IRIS100$iris012<-as.factor(IRIS100$iris012)
prop.table(table(IRIS100$iris012))
########################### cross-validation using LR ##########################
Control <- trainControl(method = "cv",number = 10)
#train function generates a candidate set of parameter values
#tuneLength controls how many are evaluated
logit.M <- train(iris012 ~., data = IRIS100, method = "glmnet", trControl = Control, tuneLength=2)
varImp(logit.M)
#Possible methods are found using 
#names(getModelInfo())

##########  classification using LR
set.seed(2307) 
indexes<-sample(nrow(IRIS100),0.8*nrow(IRIS100),replace = F) # training set 80%  100*0.8=80 entries, test set 20%  20 entries 
sort(indexes)
train<-IRIS100[indexes,]; dim(train)
test<-IRIS100[-indexes,]; dim(test)
############ LR model #######################
model=glm(iris012~Sepal.Width+Petal.Width,family=binomial,data=train)#AIC: 10
summary(model)
test_pred = predict(model, newdata = test, type = "response")
test_pred=round(test_pred,1); test_pred
write.csv(test_pred, file="test.csv")

length(test_pred)
head(test_pred)
table.LR.test=table(pred=test_pred,true=test$iris012)
table.LR.test
str(table.LR.test)
acc.LR<-round(sum(diag(table.LR.test)) / sum(table.LR.test)*100,2)
acc.LR
