#ref: https://www.jamleecute.com/hierarchical-clustering-%E9%9A%8E%E5%B1%A4%E5%BC%8F%E5%88%86%E7%BE%A4/
# https://rpubs.com/skydome20/R-Note9-Clustering
# https://rpubs.com/jiankaiwang/hclust
# https://uc-r.github.io/hc_clustering
# this code include the following unsupervised learning codes
# Hierarchical clustering, kmean, kmedoid
# optimal number of clusters: Elbow Method, Silhouette index
setwd("/Users/klng 1/R/ML")
#install.packages("cluster")
#install.packages(pkgs)
#https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/#three-popular-methods-for-determining-the-optimal-number-of-clusters
library(cluster)
IRIS=read.csv("./iris.data", header=FALSE); dim(IRIS)
head(IRIS)
# remove NA records
na.omit(IRIS); dim(IRIS) #check whether any records removed after na.omit
head(IRIS)
names(IRIS)=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species")
head(IRIS)
# unsuperivised learning --> DELETE the field 'Species' from our input data
iris=IRIS[,-5]
head(iris)
euclid = dist(x = iris, method = "euclidean")
manhatt = dist(x = iris, method = "manhattan")
# Let the graphics be presented in a layout of one line and two columns, 
# if you want to restore it, please use dev.off()

h.E.cluster <- hclust(euclid)
plot(h.E.cluster, xlab="euclidean", ylab="height",col="red")
dev.off()
h.M.cluster <- hclust(manhatt) 
plot(h.M.cluster, xlab="manhattan", ylab="height",col="blue")
########### use different clustering linkage methods #######################
## The single-linkage method has the effect of bigger dominate
# the complete linkage and average linkage has the effect of going forward together.

plot(hclust(manhatt, method="single"), xlab="manhattan", ylab="height",col="green")
plot(hclust(manhatt, method="average"), xlab="manhattan", ylab="height",col="purple")
plot(hclust(manhatt, method="complete"), xlab="manhattan", ylab="height",col="pink")
plot(hclust(manhatt, method="centroid"), xlab="manhattan", ylab="height",col="orange")
plot(hclust(manhatt, method="ward.D2"), xlab="manhattan", ylab="height",col="magenta")
abline(h=9, col="black")  # cutting the dendrogram

#par(mfrow=c(1,2))  # not working on my mac
# Agglomerative  method agnes
# agnes is very similar to hclust(), but the difference is that it can 
#additionally calculate the agglomerative coefficient
hc2 <- agnes(manhatt, method = "ward")
# Agglomerative coefficient
#The closer the aggregation coefficient to 1, the stronger the clustering structure is.
signif(hc2$ac,3)
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html
# The option "ward.D" (equivalent to the only Ward option "ward" in R versions <= 3.0.3) 
#does not implement Ward's (1963) clustering criterion, whereas option "ward.D2" implements
#that criterion (Murtagh and Legendre 2014). With the latter, the dissimilarities 
#are squared before cluster updating. 
#Note that agnes(*, method="ward") corresponds to hclust(*, "ward.D2").
# The above results can also be executed using the anges() function. 
#This function is very similar to hclust(), but the difference is that it can 
#additionally calculate the agglomerative coefficient. The aggregation coefficient 
#is a measure of the degree to which the clustering structure is identified. 
#The closer the aggregation coefficient is to 1, the stronger the clustering structure is. 
#In this case, the clustering coefficient using Euclidean distance combined with the 
#Huade connection algorithm has a performance of up to 99%.
h.cluster=hclust(manhatt, method="ward.D2")
h.cluster=hclust(euclid, method="ward.D2")
cut.h.cluster <- cutree(h.cluster, k=3)  #three clusters
cut.h.cluster   # the results of 3 clusters memberships

confusion.matrix=table(cut.h.cluster, IRIS$Species) 
confusion.matrix
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
#for the 3rd class virginica
precision.virginica= signif(confusion.matrix[3,3]/sum(confusion.matrix[3,]) , 3)
precision.virginica
sensitive.virginica= signif(confusion.matrix[3,3]/sum(confusion.matrix[,3]),3)#sum over column
sensitive.virginica

out1=c("precision.setosa","sensitive.setosa","precision.veriscolor", "sensitive.veriscolor", "precision.virginica", "sensitive.virginica")
out2=c(signif(precision.setosa,3),signif(sensitive.setosa,3),signif(precision.veriscolor,3), signif(sensitive.veriscolor,3),signif(precision.virginica,3),signif(sensitive.virginica,3))
out=rbind(out1, out2); out
# Comparison of clustering results and actual results
#It seems that this grouping is very successful in dividing setosa into 
#the first group; versicolor into the second group;
#However, virginica seems to be having a little trouble

plot(table(IRIS$Species, cut.h.cluster), main = "Confusion Matrix for Species Clustering", xlab = "Species", ylab = "Cluster")

#install.packages("ggplot"), require SAME as using library
require(factoextra)
library("ggplot2")
ggplot(data = IRIS,mapping = aes(x = Petal.Length, y = Petal.Width))+geom_point(aes(col = Species))
#manhattan distance
#[,1]               [,2]               [,3]                   [,4]                  
#out1 "precision.setosa" "sensitive.setosa" "precision.veriscolor" "sensitive.veriscolor"
#out2 "1"                "1"                "0.754"                "0.98"                
#[,5]                  [,6]                 
#out1 "precision.virginica" "sensitive.virginica"
#out2 "0.971"               "0.68"    

# euclidean distance
#[,1]               [,2]               [,3]                   [,4]                  
#out1 "precision.setosa" "sensitive.setosa" "precision.veriscolor" "sensitive.veriscolor"
#out2 "1"                "1"                "0.766"                "0.98"                
#[,5]                  [,6]                 
#out1 "precision.virginica" "sensitive.virginica"
#out2 "0.972"               "0.7"   


## using kmeans clustering RESULTS ARE NOT STALBE !!!!  ###################
kmeans.cluster=kmeans(iris,centers=3,iter.max=100,algorithm=c("Hartigan-Wong"),nstart=1)
summary(kmeans.cluster)
# variance inside the cluster
kmeans.cluster$withinss
out=table(kmeans.cluster$cluster, IRIS$Species) 
out
# acc for the three species classifier
acc=signif(sum(diag(out))/ sum(out), 3 )
acc
#for the 1st class, sentosa
#precision (also called positive predictive value)
#Precision = TP/(TP+FP) for setosa
precision.setosa= signif(out[1,1]/sum(out[1,]),3)#sum over row NOT column
precision.setosa
#recall (also known as sensitivity) is the fraction of relevant instances 
#that were retrieved #Recall (aka Sensitivity) = = TP/(TP+FN)
sensitive.setosa= signif(out[1,1]/sum(out[,1]),3)#sum over column
sensitive.setosa
#for the 2nd class veriscolor
precision.veriscolor= signif(out[2,2]/sum(out[2,]) , 3)
precision.veriscolor
sensitive.veriscolor= signif(out[2,2]/sum(out[,2]),3)#sum over column
sensitive.veriscolor
#for the 3rd class virginica
precision.virginica= signif(out[3,3]/sum(out[3,]) , 3)
precision.virginica
sensitive.virginica= signif(out[3,3]/sum(out[,3]),3)#sum over column
sensitive.virginica

library(factoextra)
library(NbClust)
require(factoextra)
fviz_cluster(kmeans.cluster, # clustering result
             data = iris, # data
             geom = c("point","text"), # point & label (point & label)
             frame.type = "norm") # Frame type

require(cluster)

# pam = Partitioning Around Medoids
kmedoid.cluster <- pam(iris, k=3) 

# Variation within the group
kmedoid.cluster$objective
# Comparison between grouping results and actual results
out.kmedoid=table(kmedoid.cluster$clustering, IRIS$Species)
out.kmedoid
# performance of kmedoid ################
# acc for the three species classifier
acc=signif(sum(diag(out.kmedoid))/ sum(out.kmedoid), 3 )
acc
#for the 1st class, sentosa
#precision (also called positive predictive value)
#Precision = TP/(TP+FP) for setosa
precision.setosa= signif(out.kmedoid[1,1]/sum(out.kmedoid[1,]),3)#sum over row NOT column
precision.setosa
#recall (also known as sensitivity) is the fraction of relevant instances 
#that were retrieved #Recall (aka Sensitivity) = = TP/(TP+FN)
sensitive.setosa= signif(out.kmedoid[1,1]/sum(out.kmedoid[,1]),3)#sum over column
sensitive.setosa
#for the 2nd class veriscolor
precision.veriscolor= signif(out.kmedoid[2,2]/sum(out.kmedoid[2,]) , 3)
precision.veriscolor
sensitive.veriscolor= signif(out.kmedoid[2,2]/sum(out.kmedoid[,2]),3)#sum over column
sensitive.veriscolor
#for the 3rd class virginica
precision.virginica= signif(out.kmedoid[3,3]/sum(out.kmedoid[3,]) , 3)
precision.virginica
sensitive.virginica= signif(out.kmedoid[3,3]/sum(out.kmedoid[,3]),3)#sum over column
sensitive.virginica
out1=c("precision.setosa","sensitive.setosa","precision.veriscolor", "sensitive.veriscolor", "precision.virginica", "sensitive.virginica")
out2=c(signif(precision.setosa,3),signif(sensitive.setosa,3),signif(precision.veriscolor,3), signif(sensitive.veriscolor,3),signif(precision.virginica,3),signif(sensitive.virginica,3))
out=rbind(out1, out2); out

# Visualize k-medoid clustering results (based on ggplot2 syntax)
#fviz_cluster = Visualize Clustering Results
require(factoextra) #library and require load and attach add-on packages.
fviz_cluster(kmedoid.cluster, # clustering result
             data = iris, # data
             geom = c("point"), # point (point)
             frame.type = "norm") # Frame type

#Elbow Method applied in hierarchical clustering
#fviz_nbclust=Dertermining and Visualizing the Optimal Number of Clusters
fviz_nbclust(iris, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) + 
  labs(title="Elbow Method for HC") +
  geom_vline(xintercept = 3, # where X=3
             linetype = 2) # draw a dotted line

# Elbow Method applied in K-Means
fviz_nbclust(iris, 
             FUNcluster = kmeans,# K-Means
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) +
  labs(title="Elbow Method for K-Means") +
  geom_vline(xintercept = 3, # where X=3
             linetype = 2) # draw a vertical dotted line

# Elbow Method applied in K-Medoid
fviz_nbclust(iris, 
             FUNcluster = pam,   # K-Medoid
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) +
  labs(title="Elbow Method for K-Medoid") +
  geom_vline(xintercept = 3, # where X=3
             linetype = 2) # draw a vertical dotted line
#When wss is about k=3, there is a turning curve, that is, the decline of wss 
#begins to decrease (tend to be stable)

#Average Silhouette Method
# Avg. Silhouette applied in K-Means
fviz_nbclust(iris, 
             FUNcluster = kmeans,   # K-Means
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +
  labs(title="Avg.Silhouette Method for K-Means") 

ss <- silhouette(kmedoid.cluster$cluster, dist(iris)); ss
mean(ss[, 3])
plot(ss)
# Gap Statistic Method
gap_stat=clusGap(x=iris,FUNcluster=hcut,nstart=25,K.max=10,B=100,spaceH0 =c("scaledPCA"))
fviz_gap_stat(gap_stat)
#According to Gap(k) statistics, the best number of clusters is 6 clusters
out=NbClust(data=iris,diss=NULL,distance="euclidean",min.nc=2, max.nc=15, method ="ward.D")
out; 
summary(out)
#Spectral Clustering
#install.packages("kernlab")
# https://www.bilibili.com/s/video/BV1nU4y1h7U8 
# http://hk.uwenku.com/question/p-xnhmndpc-sq.html
# http://www.biostatistic.net/thread-49091-1-1.html 
