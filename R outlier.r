#install.packages("neuralnet", dependencies = FALSE)
setwd("/Users/klng 1/R/ML")
iris=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
names(iris) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")

iris<- na.omit(iris) 
x1=iris[,1]; x2=iris[,2]; x3=iris[,3]; x4=iris[,4];
# outliers detection 
boxplot.stats(x1)$out 
boxplot(x1) 

df <- data.frame(x1, x2, x3, x4) 
names(df) = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
head(df) 
attach(df) 
# find the index of outliers from x 
(a <- which(x1 %in% boxplot.stats(x1)$out)) 
(b <- which(x2 %in% boxplot.stats(x2)$out)) 
(c <- which(x3 %in% boxplot.stats(x3)$out))
(d <- which(x4 %in% boxplot.stats(x4)$out))
detach(df)
# outliers in both x and y 
(outlier.list1 <- intersect(a,b)) # 0 result
(outlier.list2 <- intersect(a,c)) # 0 result
(outlier.list3 <- intersect(a,d)) # 0 result
(outlier.listb <- intersect(b,b)) # 4 result
boxplot.stats(x2)$out 
boxplot(x2) 
#plot(df) 
#points(df[outlier.list1,], col="red", pch="+", cex=2.5) 
#outliers are labeled with “+” in red


# outliers in either x or y 
(outlier.list_ab <- union(a,b)) 
plot(df) 
points(df[outlier.list_ab,], col="blue", pch="x", cex=2) 
# outliers are labeled with “x” in blue

#install.packages("DMwR")
library(DMwR) # first time not success, type one more time
library(DMwR)
# remove "Species", which is a categorical column 
iris2 <- iris[,1:4] 
outlier.scores <- lofactor(iris2, k=5) 
plot(density(outlier.scores)) 
# pick top 5 as outliers 
outliers <- order(outlier.scores, decreasing=T)[1:5] 
# who are outliers 
print(outliers) 
print(iris2[outliers,]) 
n <- nrow(iris2) 
labels <- 1:n 
labels[-outliers] <- "." 

biplot(prcomp(iris2), cex=.8, xlabs=labels) 

pch <- rep(".", n) 
pch[outliers] <- "+" 
col <- rep("black", n) 
col[outliers] <- "red" 
pairs(iris2, pch=pch, col=col) 

# remove species from the data to cluster 
iris2 <- iris[,1:4] 
kmeans.result <- kmeans(iris2, centers=3) 
# cluster centers 
kmeans.result$centers 
# cluster IDs 
kmeans.result$cluster 
# calculate distances between objects and cluster centers 
centers <- kmeans.result$centers[kmeans.result$cluster, ] 
distances <- sqrt(rowSums((iris2 - centers)^2)) 
# pick top 5 largest distances 
outliers <- order(distances, decreasing=T)[1:5] 
# who are outliers 
print(outliers) 
print(iris2[outliers,]) 



### Outlier Detection with RLOF, NOT for Microsoft OS ???
#install.packages("Rlof”)
library(Rlof)

outlier.scores <- lof(iris2, k=5) 
# try with different number of neighbors (k = 5,6,7,8,9 and 10) 
outlier.scores <- lof(iris2, k=c(5:10)) 
head(outlier.scores)
################
                 
                 

