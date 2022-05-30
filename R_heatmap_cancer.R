setwd("/Users/user/Documents/R/Hjs")
require(pheatmap)
require(RColorBrewer)

hjs3cancer <-read.csv("./Hjs3_cancer_map.csv",header=TRUE) # need to type ‘.csv’ 
#althought the filename extension does not have the ‘.csv’
dim(hjs3cancer); is.numeric(hjs3cancer) # strange, it is not numeric
M_hjs3cancer<-as.matrix(hjs3cancer)
dim(M_hjs3cancer); diag(M_hjs3cancer); hjs3cancer[0:0,] ; hjs3cancer[0:1,] 
hjs3cancer[1:1,] #SAME as above, Rstudio lists the field names as well, even we use 1:1  
hjs3cancer[1:3,] 
hjs3cancer[,0] # column ‘0’ is null
hjs3cancer[,1] # column ‘1’ are text
T_hjs3cancer<-t(hjs3cancer)
dim(T_hjs3cancer); is.numeric(T_hjs3cancer) # strange, it is numeric
#hjs3cancer<-t(T_hjs3cancer)  # transpose again, so it is numeric
is.numeric(hjs3cancer)
Sym_hjs3cancer<- hjs3cancer[1:17,1:17] + T_hjs3cancer[1:17,1:17] #row ‘0’ is text
Sym_hjs3cancer[1:2,1:2]  #graphics.off()
pheatmap(Sym_hjs3cancer, color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100), 
         fontsize_row=10, fontsize_col=10)
# below, see https://stackoverflow.com/questions/67403690/how-do-i-add-a-legend-to-my-heatmap-in-r
#pheatmap(Sym_hjs3stn, color=brewer.pal(9,"Oranges"))

hjs4cancer <-read.csv("./Hjs4_cancer_map.csv",header=TRUE) # need to type ‘.csv’ althought the filename extension does not have the ‘.csv’
dim(hjs4cancer); is.numeric(hjs4cancer) # strange, it is not numeric
hjs4cancer[16:17, 1:3] 

T_hjs4cancer<-t(hjs4cancer)
dim(T_hjs4cancer); is.numeric(T_hjs4cancer) # strange, it is numeric
#hjs4cancer<-t(T_hjs4cancer)  # transpose again, so it is numeric
is.numeric(hjs4cancer)
Sym_hjs4cancer<- hjs4cancer[1:17,1:17] + T_hjs4cancer[1:17,1:17] #row ‘0’ is text
Sym_hjs4cancer[1:2,1:2]  #graphics.off()
pheatmap(Sym_hjs4cancer, color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100), 
         fontsize_row=10, fontsize_col=10)

