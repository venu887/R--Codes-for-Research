#ref: https://www.jamleecute.com/hierarchical-clustering-%E9%9A%8E%E5%B1%A4%E5%BC%8F%E5%88%86%E7%BE%A4/
# https://rpubs.com/skydome20/R-Note9-Clustering
# https://rpubs.com/jiankaiwang/hclust
# https://uc-r.github.io/hc_clustering
#Hierarchical clustering
setwd("/Users/klng 1/R/ML")
#install.packages("cluster")
#install.packages(pkgs)
#https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/#three-popular-methods-for-determining-the-optimal-number-of-clusters
library(cluster)
setwd("/Users/klng 1/R/ML")
xyz=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
xyz=iris
dim(xyz); col=dim(xyz)[2]; col
col_minus1=col-1; col_minus1
xyz[,5]; names(xyz);
names(xyz)= c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")
head(xyz)
xyz$Species<-as.factor(xyz$Species) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(xyz$Species)),3)
fraction[1]; fraction[2]; fraction[3]
# unsuperivised learning --> delete the 5th column 'Species' from input data
XYZ=xyz[,-5]
head(XYZ); dim(XYZ)
#library(fpc)
#library(dbscan)
#ds <- dbscan(XYZ[,-1], eps=50, minPts=10) 
# compare clusters with original class labels 
#ds
#plot(ds, iris2) # there is a display problem !
#plot(ds, XYZ[c(2,4)])#

#if (is.na(XYZ)=="TRUE") {print XYZ}
XYZ=XYZ[1:100,]; dim(XYZ)
euclid = dist(x = XYZ, method = "euclidean")
manhatt = dist(x = XYZ, method = "manhattan")
# Let the graphics be presented in a layout of one line and two columns, 
# if you want to restore it, please use dev.off()

h.cluster=hclust(manhatt, method="ward.D2")
cut.h.cluster <- cutree(h.cluster, k=2)  #two clusters
cut.h.cluster   # the clustering results

confusion.matrix=table(xyz[1:100,]$Species, cut.h.cluster) 
confusion.matrix; dim(confusion.matrix)

TN=confusion.matrix[1,1]; TN
FP=confusion.matrix[1,2]; FP
FN=confusion.matrix[2,1]; FN
TP=confusion.matrix[2,2]; TP
total = TN+TP+FN+FP; total
accuracy =  (TP + TN)/total ; signif(accuracy,3)
Sensitivity  = TP/(TP+FN);signif(Sensitivity,3)
Specificity = TN/(TN+FP);signif(Specificity,3)
# performance
#for the 1st class, sentosa
#precision (also called positive predictive value)
#Precision = TP/(TP+FP) for setosa
precision.setosa= signif(confusion.matrix[1,1]/sum(confusion.matrix[1,]),3)#sum over row NOT column
precision.setosa
#recall (also known as sensitivity) is the fraction of relevant instances 
#that were retrieved #Recall (aka Sensitivity) = = TP/(TP+FN)
sensitive.setosa= signif(confusion.matrix[1,1]/sum(confusion.matrix[,1]),3)#sum over column
sensitive.setosa
#for the 2nd class veriscolor
precision.veriscolor= signif(confusion.matrix[2,2]/sum(confusion.matrix[2,]) , 3)
precision.veriscolor
sensitive.veriscolor= signif(confusion.matrix[2,2]/sum(confusion.matrix[,2]),3)#sum over column
sensitive.veriscolor


out1=c("precision.setosa","sensitive.setosa","precision.veriscolor", "sensitive.veriscolor")
out2=c(signif(precision.setosa,3),signif(sensitive.setosa,3),signif(precision.veriscolor,3), signif(sensitive.veriscolor,3))
out=rbind(out1, out2); out

# Comparison of clustering results and actual results
#It seems that this grouping is very successful in dividing setosa into 
#the first group; versicolor into the second group;
#However, virginica seems to be having a little trouble

plot(table(XYZ$Species, cut.h.cluster), main = "Confusion Matrix for Species Clustering", xlab = "Species", ylab = "Cluster")

#install.packages("ggplot")
require(factoextra)
library("ggplot2")
ggplot(data = XYZ,mapping = aes(x = Petal, y = Petal))+geom_point(aes(col = CL))

