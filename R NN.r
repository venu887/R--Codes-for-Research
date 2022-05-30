#install.packages("neuralnet", dependencies = FALSE)
setwd("/Users/klng 1/R/ML")
iris=read.csv("iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")

iris<- na.omit(iris) 
library(neuralnet)
index <- sample(2, nrow(iris), replace=T, prob=c(0.7, 0.3)) 
# training sets [70%], test set [30%], CANNOT set replace=F
nn.trainset<- iris[index==1,] 
nn.testset <- iris[index==2,] 
nn.pred<-neuralnet(Species ~ ., data=nn.trainset, hidden=10, threshold=0.01, learningrate=0.01, algorithm="backprop")
names(nn.pred)
head(nn.pred$result.matrix)
nn.test<-compute(nn.pred, nn.testset)
names(nn.test)
head(nn.test$net.result)
graphics.off()
plot(nn.pred)

nn.trainset<- iris[index==1,] 
nn.testset <- iris[index==2,] 
nn.pred<-neuralnet(Species ~ ., data=nn.trainset, hidden=10, threshold=0.01, learningrate=0.01, algorithm="backprop")
names(nn.pred)
head(nn.pred$result.matrix)
nn.test<-compute(nn.pred, nn.testset)
names(nn.test)
idx <- apply(nn.test$net.result, 1, which.max) 
predicted <- c('setosa', 'versicolor', 'virginica')[idx] 
table(actual= nn.testset$Species, pred=predicted)

