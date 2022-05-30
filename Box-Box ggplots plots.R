#BOX PLOT
rm(list = ls())
setwd("J:/V/R")
library(ggplot2)
library(reshape2)
library(lattice) 
a<-read.csv("Cancer Exp Data.csv", header = T)
names(a)
b<-log2(a[,2:23]+1)
c<-a[,9:19]
data_long<-melt(a, id="SmallRNA")
ggplot(data_long, aes(x = variable, y = value, color = SmallRNA)) +                  
  geom_boxplot()
#https://www.youtube.com/watch?v=V3Co_UFLOI4&ab_channel=StatisticsGlobe
##### Example 1
set.seed(8642)                                               # Create random data
x <- rnorm(1000)
boxplot(x)                                                   # Basic boxplot in R

##### Example 2
y <- runif(1000)                                             # Create more variables
z<-rpois(1000, 3)
data<-data.frame(values = c(x, y, z),                      # Combine variables in data frame
                   group = c(rep("x", 1000),
                             rep("y", 1000),
                             rep("z", 1000)))
head(data)                                                   # First six rows of data
boxplot(values ~ group, data)                                # Multiple boxplots in same graph
##### Example 8
data2<-data                                                # Replicate data
data2$group<-c(rep("x1", 500), rep("x2", 500),             # Modify group variable
                 rep("y1", 500), rep("y2", 500),
                 rep("z1", 500), rep("z2", 500))
boxplot(values ~ group, data2,                               # Boxplot with manual positions
        col = c("blue", "pink"),
        at = c(1, 2, 5, 6, 9, 10))
OR
ggplot(data2, aes(x = group, y = values, fill = group)) +    # Create boxplot chart in ggplot2
  geom_boxplot()

attach(a)
names(a)
boxplot(hsa.miR.1251.5p)
quantile(hsa.miR.1251.5p, probs = c(0.25, 0.5, 0.75))
boxplot(hsa.miR.1251.5p~primary_or_metastasis, main="miR exp")
boxplot(hsa.miR.1251.5p~primary_or_metastasis, 
        xlab = "primary and metastasis",
        ylab = "Normalized log counts",
        outline=FALSE,
        medcol="red",
        main = "hsa.miR.1251.5p expression",
        notch = TRUE, 
        varwidth = TRUE, 
        names = c("Primary","Metastasis")
)

set.seed(75829547)                                                  
boxplot(b[,2:11])
names(b)
data_long<-melt(b, id="primary_or_metastasis")                                             
ggplot(data_long, aes(x = variable, y = value, color = primary_or_metastasis)) +                  
  geom_boxplot()

# Install lattice package
library("lattice")                                                 
bwplot(value ~ variable, data_long)                                


#MEANS OF THE BOX PLOTS
set.seed(2967358)                                          # Create example data
data <- data.frame(data_long)

data_means <- aggregate(data$values,                       # Means by group
                        list(data$variable),
                        mean)

boxplot(data$values ~ data$group)                          # Draw boxplot in Base R
points(x = 1:nrow(data_means),                             # Add points to plot
       y = data_means$x,
       col = "red",
       pch = 16)
text(x = 1:nrow(data_means),                               # Add text to plot
     y = data_means$x - 0.15,
     labels = paste("Mean:", round(data_means$x, 1)),
     col = "red")

ggplot(data, aes(x = group, y = values)) +                 # Draw ggplot2 boxplot
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", col = "red") +  # Add points to plot
  stat_summary(fun = mean, geom = "text", col = "red",     # Add text to plot
               vjust = 1.5, aes(label = paste("Mean:", round(..y.., digits = 1))))

# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation
boxplot(primary_or_metastasis~., data=b, notch=TRUE,
        col=(c("gold","darkgreen")),
        main="Tooth Growth", xlab="Suppliment and Dose")


# https://r-graph-gallery.com/265-grouped-boxplot-with-ggplot2.html 
#Grouped boxplot with ggplot2
a<-read.csv("Cancer Exp Data.csv", header = T)
b<-a[,c(5, 9:18)]
data_long<-melt(b, id="primary_or_metastasis")
names(data_long)
# grouped boxplot
ggplot(data_long, aes(x=variable, y=value, fill=primary_or_metastasis)) + 
  geom_boxplot()
# One box between Primary vs Metastasis 
ggplot(data_long, aes(x=variable, y=value, fill=primary_or_metastasis)) + 
  geom_boxplot() +
  facet_wrap(~primary_or_metastasis)
# one box per Primary vs Metastasis
ggplot(data_long, aes(x=variable, y=value, fill=SmallRNA)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free")



###Multiple BAR plots
# creating multiple bar plots in R
# creating plot using the above data
#Rotate X axis labels in ggplot with 90 Degree Angle
p<-ggplot(data_long,aes(x=variable, y=value, fill=SmallRNA)) + 
  geom_bar(stat="identity", position = "dodge") +
  scale_y_continuous("Normalized TMM Values", expand = c(0, 0)) +
  scale_x_discrete("miRNA")+
  labs(title="OV-Three Cell lines Expression values(Log(TMM+1))")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, vjust = 0),
        axis.line = element_blank(),
        axis.ticks.x = element_blank())+
  coord_flip() # This addition of code indicates flip our original figure to horizontal 
print(p)          
