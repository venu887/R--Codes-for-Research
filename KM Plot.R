# Survival model KM
#install.packages("survminer")
setwd("J:/V/R")
rm(list = ls())
library(limma)
library(glmnet)
library(survival)
library(ggplot2) 
library(dplyr)
#Surv(): Creates a survival object.
#survfit(): Fits a survival curve using either a formula, of from a previously fitted Cox model.
#coxph(): Fits a Cox proportional hazards regression model.
#cox.zph(): Tests the proportional hazards assumption of a Cox regression model.
#survdiff(): Tests for differences in survival between two groups using a log-rank / Mantel-Haenszel test
rm(list = ls())
# Convert raw reads to Z scores

a<-read.csv("TCGA3GEOCU_match.csv", row.names = 1)
a<-read.csv("Cancer RSEM.csv", row.names = 1)
y <- Surv(a$os, a$event)
km.model<-survfit(Surv(dat$os,dat$event) ~dat$hsa.mir.526b.1, data = dat, type="kaplan-meier")
print(km.model)
plot(km.model)
summary(km.model)
summary(km.model,times=c(25,50,75,100,125))

library(survminer)
ggsurvplot(km.model, data = a, risk.table = T, conf.int = T, ggtheme = theme_minimal())

Km.model2<-survfit(Surv(a$os, a$event)~stage, type="kaplan-meier", data = a)
print(Km.model2)
summary(Km.model2)
summary(km.model,times=c(25,50,75,100,125))
ggsurvplot(Km.model2, data=a, pval=T, risk.table=T,conf.int = T, pval.method = T, 
           ggtheme =theme_minimal() )
#estimate COX PH-model-cox regression models
cox_reg1<-coxph(Surv(a$os, a$event)~stage, data = a)
summary(cox_reg1)

# age must be yes or no with medium 
cox_reg2<-coxph(Surv(a$os, a$event)~ hsa.mir.137+hsa.mir.2116+hsa.mir.1227+hsa.mir.548f.1+hsa.mir.1251+hsa.mir.1225, data = a)
summary(cox_reg2)
#we can create formula from summary results 
# insert expression of miR in formula
hsa.mir.137=2
hsa.mir.548f.1=1.5
hsa.mir.1225=2.5
log_overall_risk=1.2438*hsa.mir.137 +0.20153*hsa.mir.548f.1-0.44265*hsa.mir.1225
print(log_overall_risk)
exp(log_overall_risk)
#grand-mean center continuous covariates
a$hsa.mir.137<-scale(a$hsa.mir.137, center = T, scale = F)
a$hsa.mir.1225<-scale(a$hsa.mir.1225, center = T, scale = F)
a$hsa.mir.548f.1<-scale(a$hsa.mir.548f.1, center = T, scale = F)
#change reference stages re-ordering
a$stage<-factor(a$stage,levels = c("es", "ls"))

#estimate cox proportional hazard model with  categorical and continuous covariates
cox_reg3<-coxph(Surv(a$os, a$event)~stage+ hsa.mir.137+hsa.mir.548f.1+hsa.mir.1225, data = a)
summary(cox_reg3)

#nested comparison 
library(tidyr)
##estimate cox proportional hazard model- COX PH-model-cox regression models
cox_reg1<-coxph(Surv(a$os, a$event)~stage, 
                data = drop_na(a,os,event,hsa.mir.137,hsa.mir.2116,hsa.mir.1227,hsa.mir.548f.1,hsa.mir.1251,hsa.mir.1225))
##estimate cox proportional hazard model- COX PH-model-cox regression models
cox_reg2<-coxph(Surv(a$os, a$event)~ stage + hsa.mir.137+hsa.mir.548f.1+hsa.mir.1225, 
                data = drop_na(a,os,event,hsa.mir.137,hsa.mir.2116,hsa.mir.1227,hsa.mir.548f.1,hsa.mir.1251,hsa.mir.1225))
anova(cox_reg1, cox_reg2) 


#confidencial Interval True is good one to know
plot(km.model, conf.int = F, xlab = "Time", ylab = "% Alive=S(t)", main="KM-Model", las=1)
plot(km.model, conf.int = T, xlab = "Time", ylab = "% Alive=S(t)", main="KM-Model", las=1)
abline(h=0.5, col="red")
plot(km.model, 
     conf.int = T, xlab = "Time", 
     ylab = "% Alive=S(t)", 
     main="KM-Model", 
     las=1, mark.time = TRUE)

Km.Model2<-survfit(y~a$FAM83D)
plot(Km.Model2, conf.int = T, xlab = "Time", ylab = "% Alive=S(t)", main="KM-Model", las=1)

d<-read.csv("TCGA3GEOCU_match.csv")
fit<- survfit(Surv(d$os,d$event) ~d$ACSL1, data = d)
plot(fit)
e<-ggsurvplot(fit, data = d, pval = T,
           legend.labs = c("Low Exp ", "High exp"), # divide expression patients high vs low
           xlab = "time(months)",
           ylab = "surviVal time",
           legend.title = "Expression") 


