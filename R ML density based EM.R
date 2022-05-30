################ Density-based Clustering ##############
setwd("/Users/klng 1/R/ML")
xyz=read.csv("./iris.data",header=FALSE) #CANNOT use header=TRUE, otherwise
dim(xyz); col=dim(xyz)[2]; col
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
#install.packages("fpc")
#install.packages("dbscan")
library(fpc)
library(dbscan)
ds <- dbscan(XYZ, eps=0.42, minPts=5) 
# number of minimum points in the eps region (for core points). 
# Default is 5 points.
# compare clusters with original class labels 
ds
plot(ds, XYZ[c(1,4)]) # extract rows 1 and 4 data from XYZ 
plotcluster(XYZ, ds$cluster) 
#the data are projected onto each classes

# compare clusters with original class labels
confusion.matrix=table(xyz$Species,ds$cluster)
confusion.matrix
#In the above data table, 1 to 3 are the three identified clusters, 
#and 0 means noisy data, that is, 
#objects that do not belong to any cluster

# how to identify optimal eps values ? ###################
d=dist(XYZ)  # XYZ consists of values of the 4 measurements
library(ggplot2)
interval=cut_interval(d,50) #makes n=50 groups with equal range
table(interval)
# (0.709,0.85] interval has the highest frequency of occurrence
#set eps between 0.7 to 0.9

EPS =0;  min = 0
for (EPS in 7:9) {
  for ( min in 5:6 ) {
     ds <- dbscan(XYZ, eps=EPS/10, minPts=min)
     print (ds)
     plot(ds, XYZ[c(1,4)]); print (EPS); print(min)
                      }  
                  }


##   using EM algorithm ##########
library(mclust)
em = Mclust(XYZ); em
summary(em, parameters=TRUE)
plot(em)
em.BIC=mclustBIC(XYZ)
em.BIC # top three models based on BIC
em.summary=summary(em.BIC, data=XYZ)
em.summary
em
plot(em$data)
plot(em$BIC);
plot(em$classification)
plot.Mclust(em,what = "density")
plot.Mclust(em,what = "classification")
mclust2Dplot(XYZ[,1:2], classification = em.summary$classification, parameters = em.summary$parameters, col="red")
mclust2Dplot(XYZ[,3:4], classification = em.summary$classification, parameters = em.summary$parameters, col="black")

em.density=densityMclust(XYZ[,1:2])
plot(em.density, XYZ[,1:2], col="grey",nlevels=55)
plot(em.density, XYZ[,1:2], type="presp", col=gery(0.8))  # not working

em.density34=densityMclust(XYZ[,3:4])
plot(em.density34, XYZ[,3:4], col="grey",nlevels=55)
plot(em.density34, XYZ[,3:4], type="presp", col=gery(0.8))  # not working
