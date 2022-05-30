# using unsupervised learning algorithms to classify cancer cohort
# clValid: an R package for cluster validation
setwd("/Users/klng 1/R/ML")
library(cluster)
require(factoextra)
library("ggplot2")
xyz=read.csv("./abc_demo.csv",header=TRUE)
dim(xyz); col=dim(xyz)[2]; col
col_minus1=col-1; col_minus1
xyz[,7]; names(xyz);
xyz$CL
xyz$CL<-as.factor(xyz$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(xyz$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
#rows 2:23 [22 samples] stages1+2, rows 24-301 stages3+4
# data imbalance problem # size of '0' >> size of '1'
dim(xyz)
xyz2=na.omit(xyz); dim(xyz2)

# unsuperivised learning --> delete the 7th column 'CL' from input data
xyz3=xyz2[,-7]; head(xyz3)
XYZ=xyz3[,-1]
head(XYZ); dim(XYZ)

euclid = dist(x = XYZ, method = "euclidean")
#Ward method forms clusters by maximizing within- clusters homogeneity
h.E.cluster=hclust(euclid,method="ward.D2")
# method="centroid"
plot(h.E.cluster, xlab="Euclidean", ylab="height",col="red")
cut.h.cluster <- cutree(h.E.cluster, k=2)  # TWO clusters
cut.h.cluster   # the results of 2 clusters memberships
dev.off()

hc2 <- agnes(euclid, method = "ward")
# Agglomerative coefficient
#The closer the aggregation coefficient to 1, the stronger the clustering structure is.
signif(hc2$ac,3)

#confusion.matrix=table(cut.h.cluster, xyz$CL) 
#confusion.matrix

#Elbow Method applied in hierarchical clustering
#fviz_nbclust=Dertermining and Visualizing the Optimal Number of Clusters
fviz_nbclust(XYZ, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
k.max = 12 )         # max number of clusters to consider
labs(title="Elbow Method for HC") 
#  geom_vline(xintercept = 3, # where X=3
#             linetype = 2) # draw a dotted line


# pam = Partitioning Around Medoids
kmedoid.cluster <- pam(XYZ, k=2, metric = c("euclidean")) 
# Variation within the group
kmedoid.cluster$objective
# Comparison between grouping results and actual results
table.testset=table(true=xyz$CL,pred=kmedoid.cluster$clustering)
table.testset

TN=table.testset[1,1]; TN
TP=table.testset[2,2]; TP
FN=table.testset[2,1]; FN
FP=table.testset[1,2]; FP
#accuracy, Q, sensitivity, SN, specificity, SP, F1 measure, F1, positive predictive rate, 
#STP , negative predictive rate, STN , false positive rate (alpha), and 
#false negative rate (beta), they are defined as 
Q = (TP+TN)/(TP+TN+FP+FN); signif(Q,3)
SN = TP/(TP+FN); signif(SN,3)
SP = TN/(TN+FP); signif(SP,3)
STP = TP/(TP+FP); signif(STP,3)
F1 =1/( 1/2*(1/SN + 1/STP) ); signif(F1,3)
#negative predictive rate,
STN = TN/(TN+FN); signif(STN,3)
MCC = ((TP+TN)-(FP*FN))/ ( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ); signif(MCC,3)
#False positive rate (alpha),
alpha= FP / (FP + TN) #= 1 – specificity
signif(alpha,3)
#False negative rate (beta),
beta= FN / (TP + FN); signif(beta,3)


out1=c("Q","SN","SP", "STP", "F1", "STN","MCC","alpha","beta")
out2=c(signif(Q,3),signif(SN,3),signif(SP,3), signif(STP,3),signif(F1,3),signif(STN,3),signif(MCC,3),signif(alpha,3),signif(beta,3))
out=rbind(out1, out2); out

#Average Silhouette Method
# Avg. Silhouette applied in kmedoid
fviz_nbclust(XYZ, 
             FUNcluster = pam,   # 
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +
  labs(title="Avg.Silhouette Method for kmedoid") 
ss <- silhouette(kmedoid.cluster$cluster, dist(XYZ)); ss
mean(ss[, 3])
plot(ss)
clusplot(XYZ, kmedoid.cluster$clustering)

#Average Silhouette Method
# Avg. Silhouette applied in hclust
fviz_nbclust(XYZ, 
             FUNcluster = hcut,   # 
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +
  labs(title="Avg.Silhouette Method for hclust") 

ss <- silhouette(cut.h.cluster, dist(XYZ)); ss
mean(ss[, 3])
plot(ss)

# heatmap
x=as.matrix(xyz[, 2:6])
row.color <- rainbow(nrow(x), start=0, end=.5) 
col.color <- rainbow(ncol(x), start=0, end=.5)
hv <- heatmap(x, col = cm.colors(256), scale="column",
        RowSideColors = row.color, ColSideColors = col.color, margins=c(5,10),
        main = "Heatmap of patients")
#xlab = "DEGs", ylab= "patients", main = "Heatmap of patients")


# Gap Statistic for Estimating the Number of Clusters
#gskmn <- clusGap(XYZ, FUN = pam, nstart = 20, K.max = 8, B = 60)
