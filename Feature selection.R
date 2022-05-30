#install.packages("caret") #install.packages("caret", dependencies = c("Depends", "Suggests"))#answer 'no', if you answer 'Yes' to compiltion [take a long time]
rm(list = ls())
setwd("J:/V/R")
#Feature selection
a<-read.csv("cesc miR.csv", row.names = 1)
b<-read.csv("DEGs DEGseq2 cesc miR.csv", row.names = 1)
names(b)
rownames(a)
c<-a[str_detect(rownames(a), "hsa-mir-514b|hsa-mir-144|hsa-mir-514-2|hsa-mir-9-3|hsa-mir-451|hsa-mir-506|hsa-mir-514-3|hsa-mir-514-1|hsa-mir-509-2|hsa-mir-509-1|hsa-mir-486|hsa-mir-508|hsa-mir-509-3"),]
d<-a[str_detect(rownames(a),"hsa-mir-449b|hsa-mir-449a|hsa-mir-526b|hsa-mir-516a-1|hsa-mir-522|hsa-mir-935|hsa-mir-934|hsa-mir-519a-2|hsa-mir-516a-2|hsa-mir-1298|hsa-mir-1911|hsa-mir-215|hsa-mir-10a|hsa-mir-204|hsa-mir-519a-1|hsa-mir-375|hsa-mir-192|hsa-mir-194-1|hsa-mir-194-2"),]
e<-rbind(c,d)
ct<-t(e)

write.csv(ct,"cesc miR ML.csv")

UPdeg<-b[b$log2FoldChange>= 1 & b$padj<0.05,] #UPREGULATED GENES
Downdeg<-b[b$log2FoldChange<= -1 & b$padj<0.05,] #DOWNREGULATED GENES
b_final<-rbind(UPdeg,Downdeg)
rownames(b_final)
library(stringr)

b<-merge(OV_20,b_final, by = "row.names", all.x = FALSE,all.y = FALSE)
b[,1] # after merrging we can see row names is in first column 
rownames(b)<-b[,1] # we have to replace row names 
b[,1]
View(b)
b<-b[,-1]
new_b<-b[,1:288]
bt<-t(new_b)
write.csv(rownames(bt), "Class Label.csv")
# arrange in excel change patient names as row names 
# add new column Class Label m0 as negative "0" m1 as Positive "1"
CL<-read.csv("Class Label.csv", row.names = 1)
View(CL)
class(CL)
d<-cbind(bt,CL)
colnames(d)
write.csv(d, "OV for svm.csv")

