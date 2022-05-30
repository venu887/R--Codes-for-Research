#kernel Multiple linear regression model 
setwd("J:/V/R")
ov<-read.csv("OV for svm.csv")
pairs(ov[2:7])
d<-ov
d[,1]
rownames(d)<-d[,1]
d<-d[,-1]
dt<-t(d)
rownames(dt)
model1<-lm(d~CL+HOXA5.3202)
class(dt)

set.seed(2323)
indexes<-sample(nrow(d),0.7*nrow(d),replace = F)
indexes
sort(indexes)
train<-d[indexes,]; dim(train)
class(train)
new_train<-as.data.frame(train)
test<-d[-indexes,]; dim(test)
class(test)
test<-as.data.frame(test)
d$CL = as.factor(d$CL)
model=glm(CL~.,family=binomial,data=train)#AIC: 10
summary(model)
model_pred=zapsmall( predict(model, newdata = new_train, type = "response") )
View(model_pred)
model_pred

test_pred = predict(model, newdata = test, type = "response")
test_pred=round(test_pred,1); test_pred
View(test_pred)
length(test_pred)


