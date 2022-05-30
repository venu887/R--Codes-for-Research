# using unsupervised learning algorithms to classify cancer cohort
# clValid: an R package for cluster validation
rm(list = ls())
setwd("J:/V/R")
library(cluster)
require(factoextra)
library("ggplot2")
V=read.csv("Cancer RSEM.csv")
dim(V); col=dim(V)[2]; col
col_minus1=col-1; col_minus1
V[,17]; names(V);
V$CL
V$CL<-as.factor(V$CL) 
# IMPORTANT to change it to factor, else problematic later !!!
fraction=signif(prop.table(table(V$CL)),3)
fraction[1]; fraction[2]
# data imbalance problem # size of '0' >> size of '1'
#rows 2:23 [22 samples] stages1+2, rows 24-301 stages3+4
# data imbalance problem # size of '0' >> size of '1'
dim(V)
V2=na.omit(V); dim(V2)

# unsuperivised learning --> delete the 7th column 'CL' from input data
V3=V2[,-17]; head(V3)
V=V3[,-1]
head(V); dim(V)

euclid = dist(x = V, method = "euclidean")
#Ward method forms clusters by maximizing within- clusters homogeneity
h.E.cluster=hclust(euclid,method="ward.D2")
# method="centroid"
plot(h.E.cluster, xlab="Euclidean", ylab="height",col="red")
cut.h.cluster <- cutree(h.E.cluster, k=2)  # TWO clusters
cut.h.cluster   # the results of 2 clusters memberships
#dev.off()

hc2 <- agnes(euclid, method = "ward")
# Agglomerative coefficient
#The closer the aggregation coefficient to 1, the stronger the clustering structure is.
signif(hc2$ac,3)

#confusion.matrix=table(cut.h.cluster, V$CL) 
#confusion.matrix

#Elbow Method applied in hierarchical clustering
#fviz_nbclust=Dertermining and Visualizing the Optimal Number of Clusters
fviz_nbclust(V, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 12 )         # max number of clusters to consider
labs(title="Elbow Method for HC") 
#  geom_vline(xintercept = 3, # where X=3
#             linetype = 2) # draw a dotted line


# pam = Partitioning Around Medoids
kmedoid.cluster <- pam(V, k=2, metric = c("euclidean")) 
# Variation within the group
kmedoid.cluster$objective
# Comparison between grouping results and actual results
table.testset=table(true=V2$CL,pred=kmedoid.cluster$clustering)
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
alpha= FP / (FP + TN) #= 1 - specificity
signif(alpha,3)
#False negative rate (beta),
beta= FN / (TP + FN); signif(beta,3)


out1=c("Q","SN","SP", "STP", "F1", "STN","MCC","alpha","beta")
out2=c(signif(Q,3),signif(SN,3),signif(SP,3), signif(STP,3),signif(F1,3),signif(STN,3),signif(MCC,3),signif(alpha,3),signif(beta,3))
out=rbind(out1, out2); out

#Average Silhouette Method
# Avg. Silhouette applied in kmedoid
fviz_nbclust(V, 
             FUNcluster = pam,   # 
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +
  labs(title="Avg.Silhouette Method for kmedoid") 
ss <- silhouette(kmedoid.cluster$cluster, dist(V)); ss
mean(ss[, 3])
plot(ss)
clusplot(V, kmedoid.cluster$clustering)

#Average Silhouette Method
# Avg. Silhouette applied in hclust
fviz_nbclust(V, 
             FUNcluster = hcut,   # 
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +
  labs(title="Avg.Silhouette Method for hclust") 

ss <- silhouette(cut.h.cluster, dist(V)); ss
mean(ss[, 3])
plot(ss)

# heatmap
x=as.matrix(V[, 2:15])
row.color <- rainbow(nrow(x), start=0, end=.5) 
col.color <- rainbow(ncol(x), start=0, end=.5)
hv <- heatmap(x, col = cm.colors(256), scale="column",
              RowSideColors = row.color, ColSideColors = col.color, margins=c(5,10),
              main = "Heatmap of patients")
#xlab = "DEGs", ylab= "patients", main = "Heatmap of patients")


# Gap Statistic for Estimating the Number of Clusters
#gskmn <- clusGap(V, FUN = pam, nstart = 20, K.max = 8, B = 60)
