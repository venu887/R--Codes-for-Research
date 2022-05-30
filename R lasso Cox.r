# Scripts are used in prognostic and diagnostic model development 
#library(edgeR)
rm(list = ls())
library(limma)
library(glmnet)
library(survival)
#cohort="BLCA" ###txt --> csv, delete 2nd row
#cohort="BRCA"  
setwd("J:/V/R")
#setwd(paste0("/Users/klng 1/R stuff/pancancer/",cohort))#
#dat<-read.csv("TCGA3GEOCU_match.csv", row.names = 1)
dat<-read.csv("Cancer RSEM.csv", row.names = 1)
any(is.na(dat)) # must be False
set.seed(1234)
# Prognostic gene with Elastic net algorithm selection 
y <- Surv(dat$os, dat$event) # os-overall survival, event-overall survival status
y <- Surv(dat$OS, dat$Event) # os-overall survival, event-overall survival status
sum(y<=0)
View(y)
x<-model.matrix(y~., dat[,c(-1:-2)])
fit <- glmnet(x,y, family = "cox", alpha = 0.5)
plot(fit)
cv.fit <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 0.5) # do not use na.omit here
plot(cv.fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
pp=coef(cv.fit, s = "lambda.min") # extracting selecting prognostic genes
pp

install.packages("msaenet")
library(msaenet)
help(msaenet)

# Prognostic gene with Lasso algorithm selection 
y <- Surv(dat$os, dat$event) 
x<-model.matrix(y~., dat[,c(-1:-2)])
cv.fit <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 1)
plot(cv.fit)
fit <- glmnet(x,y, family = "cox", alpha = 1)
plot(fit)
cv.fit$lambda.min
Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
TT=coef(cv.fit, s = "lambda.min") # extracting selecting prognostic genes
TT
#addaptive lasso
cv.fit <- cv.glmnet(x,y, family="cox", nfold = 10, alpha= 0) # rigid Cox Model
plot(cv.fit)
fit <- glmnet(x,y, family = "cox",
              nfold = 10, alpha= 0)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
aa=coef(cv.fit, s = "lambda.min")
aa
# addaptive lasso
## Extract coefficients at the error-minimizing lambda
cv.fit$lambda.min
coef(cv.fit, s = "lambda.min")
## The intercept estimate should be dropped.
best_ridge_coef <- as.numeric(coef(cv.fit, s = "lambda.min"))
best_ridge_coef
## Perform adaptive LASSO
alasso1 <- glmnet(x,y, family = "cox", alpha = 1,
                  penalty.factor = 1 / abs(best_ridge_coef))
plot(alasso1, xvar = "lambda")

## Perform adaptive LASSO with 10-fold CV
alasso1_cv <- cv.glmnet(x,y, family = "cox",
                        ## type.measure: loss to use for cross-validation.
                        ## K = 10 is the default.
                        nfold = 10,
                        ## 'alpha = 1' is the lasso penalty, and 'alpha = 0' the ridge penalty.
                        alpha = 1,
                        penalty.factor = 1 / abs(best_ridge_coef),
                        ## prevalidated array is returned
                        keep = TRUE)
## Penalty vs CV MSE plot
plot(alasso1_cv)

## Extract coefficients at the error-minimizing lambda
alasso1_cv$lambda.min
## s: Value(s) of the penalty parameter 'lambda' at which
##    predictions are required. Default is the entire sequence used
##    to create the model.
best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
best_alasso_coef1 # extracting prognostic gene selected by Adaptive lasso 
# Best subset Cox regression model construction 
# Combine all gene selected by three algorithms and use glumti package to identify which combination of genes would produce optimal prognostic signture 



#best subset regression analysis 
library(survival)
#install.packages("glmulti")
library(glmulti)
dat <- within(dat, {
  survival.vector    <- Surv(dat$os, dat$event)})
#*************************************************************
C12orf36.283422+C1S.716+CCL19.6363+COL22A1.169044+CPZ.8532+CTLA4.1493+GFPT2.9945+IQGAP3.128239+KIF20A.10112+LOC84856.84856+MIR155HG.114614+NCAPG.64151+NEK2.4751+NUF2.83540+RUFY4.285180+TROAP.10024

names(dat)
CESC<-hsa.mir.522+hsa.mir.526b+hsa.mir.935+hsa.mir.144+hsa.mir.192+hsa.mir.204+hsa.mir.508+hsa.mir.486
glmulti.coxph.out <-glmulti(survival.vector ~C12orf36.283422+C1S.716+CCL19.6363+COL22A1.169044+CPZ.8532+CTLA4.1493+GFPT2.9945+IQGAP3.128239+KIF20A.10112+LOC84856.84856+MIR155HG.114614+NCAPG.64151+NEK2.4751+NUF2.83540+RUFY4.285180+TROAP.10024, data = dat,
                            level = 1,               # No interaction considered
                            method = "h",            # Exhaustive approach
                            crit = "bic",            # AIC as criteria
                            confsetsize =10,         # Keep 5 best models
                            plotty = F, report = F,  # No plot or interim reports
                            fitfunction = "coxph")   # coxph function

#glmulti.coxph.out <-glmulti(survival.vector ~FAM83D +	CDC20 +	TPX2 + 	LECT2 + ANXA10 +	DNASE1L3 +	PON1 +	CD5L +	CYP2C9 +	ADH4 +	CFHR3 +	GHR + LCAT, data = dat,
                            level = 1,               # No interaction considered
                            method = "h",            # Exhaustive approach
                            crit = "aic",            # AIC as criteria
                            confsetsize = 5,         # Keep 5 best models
                            plotty = F, report = F,  # No plot or interim reports
                            fitfunction = "coxph")   # coxph function

#glmulti.coxph.out <-glmulti(survival.vector ~hsa.mir.2113+hsa.mir.891b+hsa.mir.1225+hsa.mir.122+hsa.mir.137+hsa.mir.2117+hsa.mir.2115+hsa.mir.376a.2+hsa.mir.490+hsa.mir.1197+hsa.mir.599+hsa.mir.323+hsa.mir.1224+hsa.mir.138.1+hsa.mir.2116+hsa.mir.483+hsa.mir.323b+hsa.mir.1227+hsa.mir.499+hsa.mir.548f.1+hsa.mir.1228+hsa.mir.548y+hsa.mir.889+hsa.mir.211+hsa.mir.2110+hsa.mir.184+hsa.mir.2114+hsa.mir.1226+hsa.mir.219.2+hsa.mir.1229+hsa.mir.675+hsa.mir.1251, data = dat,
                            level = 1,               # No interaction considered
                            method = "h",            # Exhaustive approach
                            crit = "aic",            # AIC as criteria
                            confsetsize = 5,         # Keep 5 best models
                            plotty = F, report = F,  # No plot or interim reports
                            fitfunction = "coxph")   # coxph function
## Show result for the best model
summary(glmulti.coxph.out@objects[[1]])
## Show 5 best models (Use @ instead of $ for an S4 object)
glmulti.coxph.out@formulas
res=glmulti.coxph.out
res
plot(res, labels(T))
coef(res)
plot(res, type="s")
library(ggplot2)
qqPlot(res, main="QQ Plot")
summary(res@objects[[1]])
print(res)
top <- weightable(res)
top
# univariant Cox regression
res.cox <- coxph(Surv(dat$os, dat$event) ~ dat$hsa.mir.137, data = dat)
res.cox
hsa.mir.137+hsa.mir.2116+hsa.mir.1227+hsa.mir.548f.1+hsa.mir.1251+hsa.mir.1225
roc(Surv(dat$os, dat$event) ~ hsa.mir.1225, data = dat)
#ROC CURVES
#install.packages("survivalROC")
library(knitr)
library(survivalROC)
library(pROC)
cutoff <- 125
set.seed(420)
z<-survivalROC(dat$os, dat$event, dat$FAM83D, entry = NULL, predict.time=cutoff, 
               cut.values =NULL, method = "KM", lambda = NULL, span = NULL, window ="symmetric")
plot(z$FP, z$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(z$AUC,3)), 
     ylab="TP",main="FAM83D, Method = KM \n Year = 1")
abline(0,1)

library(pROC)
ROC <- roc(dat$os~dat$event,plot=TRUE,print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")
gene1_ROC <- roc(dat$os ~ dat$FAM83D,plot=TRUE,print.auc=TRUE,col="blue",lwd = 4,print.auc.y=0.4,legacy.axes=TRUE,add = TRUE)
gene2_ROC <- roc(dat$os ~ dat$TPX2,plot=TRUE,print.auc=TRUE,col="red",lwd = 4,print.auc.y=0.4,legacy.axes=TRUE,add = TRUE)
gene3_ROC <- roc(dat$os ~ dat$ADH4,plot=TRUE,print.auc=TRUE,col="black",lwd = 4,print.auc.y=0.4,legacy.axes=F,add = TRUE)
gene4_ROC <- roc(dat$os ~ dat$LECT2,plot=TRUE,print.auc=TRUE,col="gray",lwd = 4,print.auc.y=0.4,legacy.axes=F,add = TRUE)

Genes_ROC<-roc(dat$CYP4A22, dat$CYP4A11, dat$APOA5)
print()
legend("bottomright",legend=c("os_ROC","gene1_ROC"),col=c("green","blue"),lwd=4)

GLM <- glm(dat$event~FAM83D+CDC20, family=binomial(logit), data=dat)
roc.data <-roc(dat$event~FAM83D+CDC20,data=dat)
oc.test (roc.data$FAM83D, roc.data$CDC2)
