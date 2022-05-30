#create a function 
print(seq(32,44))
print(mean(25:82))
print(sum(41:68))
# function with an argument "a"
new_function<-function(a){
  for(i in 1:a){
    b<-i^2
    print(b)
  }
}
#function without an argument 
new.function<-function(){
  for(i in 1:5){
    print(i^2)
  }
}
new.function()

#function with multiple arguments
new.function<-function(a,b,c){
  results<-a*b+c 
  print(results)
}
new.function(5,11,3)
new.function(a = 11, b = 5, c = 3)
#Calling a Function with Default Argument
new.function<-function(a=3, b=6){
  result<-a*b
  print(result)
}
new.function()
new.function(9,5)

# Lazy Evaluation of Function Create a function with arguments.
new.function <- function(a, b) {
  print(a^2)
  print(a)
  print(b)
}

# Evaluate the function without supplying one of the arguments.
new.function(6)

#Dealing with Outliers 
remove_outliers <- function(dat1, na.rm = TRUE) {
  qnt <- quantile(dat1, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(dat1, na.rm = na.rm)
  y <- dat1
  y[dat1 < (qnt[1] - H)] <- NA
  y[dat1 > (qnt[2] + H)] <- NA
  print(y)
}
y <- remove_outliers(dat1$hsa.mir.526b.1,)
boxplot(dat1, outline = F)
boxplot(y)

# REMOVE OUTLIERS AND PERFORM COX REGRESSION ANALYSIS
#Validation of miRNA prognostic power in hepatocellular carcinoma using expression data of independent datasets
# We determined each percentile of miRNA expression between the lower and upper quartiles of expression as a cutof point to divide 
#patients into high and low expression groups as described previously45. Because of the low sample number, 
#a cutof outside the lower or upper quartile of expression could result in unreliable results. 
#Afer this, the Cox regression analysis was performed separately for each cutof. 
library(survival)
setwd("J:/V/R")
rm(list = ls())
dat1<-read.csv("Cancer RSEM.1.csv", header = 1)
str(dat1)
summary(dat1)
# Required PACKAGES
library(beeswarm)
library(ggplot2)
library(stringr)
library(proxyC)
# ==== SURVIVAL FUNCTION ====
bestcutoff <- function(datavector, clintable) {
  breaks <- quantile(datavector, probs = seq(0.25, 0.75, by= 0.01))
  cutoff.table <- t(sapply(breaks, function(z) cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  #cutoff.table
  cutoff.table[order(cutoff.table[, 2]), "cutoff"][1]
}

#removing outlines in single miR
dat1<-read.csv("Cancer RSEM.csv", header = 1)
z<-dat1[colSums(dat1[,4:1046]==0)<=20,]
y<-dat1[str_detect(names(dat1), "OS|Event|hsa.mir.526b")] 
summary(dat1)
c<-dat1[dat1$hsa.mir.526b,]
#Outliers formula-1
z <- c[c$hsa.mir.526b > quantile(c$hsa.mir.526b, .10) - 1.5*IQR(c$hsa.mir.526b) & 
            c$hsa.mir.526b< quantile(c$hsa.mir.526b, .90) + 1.5*IQR(c$hsa.mir.526b), ] #rows
#Outliers formula-1
newdata <- subset(dat1,!(dat1$hsa.mir.526b > quantile(dat1$hsa.mir.526b, probs=c(.01, .99))[2] | dat1$hsa.mir.526b < quantile(dat1$hsa.mir.526b, probs=c(.01, .99))[1]) ) 
#Outliers formula-1
#outliers<-function(x){
  iqr<-IQR(x)
  q1<-as.numeric(quantile(x,0.25))
  q3<-as.numeric(quantile(x,0.75))
  mild_low<-q1-(1.5*iqr)
  mild_high<-q3+(1.5*iqr)
  new_variable<-x[x>mild_low & x<mild_high]
  return(new_variable)
  print(new_variable)
}

dat2<-newdata[,-1]
dat2<-dat1[,-1]
y<-dat2$OS
y<-as.numeric(y)
x<-dat2[,-1:-2]
event<-dat2$Event
event<-as.numeric(event)
ans = apply(x,2,function(x,y,event)coxph(Surv(y,event)~x),y=y,event=event)

# Extract data 
univ_results <- lapply(ans,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=4)
                         wald.test<-signif(x$wald["test"], digits=4)
                         beta<-signif(x$coef[1], digits=4);#coeficient beta
                         HR <-signif(x$coef[2], digits=4);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
result<-as.data.frame(res)
write.csv(result,"CESC cox uni.csv")
#write.csv(result,"OV-univariate Cox Regression Analysis.csv")




