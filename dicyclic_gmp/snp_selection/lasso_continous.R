rm(list=ls(all=TRUE))

library(glmnet)
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp_data.R')
tg<-t(snp.by.gene.uniq[[49]])#2 dim 3 by 30
t<-m.cdgmp
x<-tg#for testing purposes
#c.lasso<-function(x){
  split<-sample(nrow(x), floor(0.7*nrow(x)))
  train.d<-as.matrix(cbind(x,t)[split,])
  test.d<-as.matrix(cbind(x,t)[-split,])
  train.m<-train.d[,-ncol(train.d)]
  target<-train.d[,ncol(train.d)]
  dim(train.m)
  dim(train.d)
  length(target)
  target<-abs(target)
  fit <- glmnet(train.m, target, family = "poisson")
  target
  pred
  ?glmnet
  fit$beta
  plot(fit$beta)
  plot(fit)
  
  plot(fit$a0)
  #summary(fit)
  
  fit$beta
  cv<-cv.glmnet(train.m,target,nfolds = 10)#cross_validate to get the best lambda
  plot(cv)
  plot(fit, xvar = "lambda", label = TRUE)
  plot(fit, xvar = "dev", label = TRUE)
  
  pred<-predict(fit, train.m, type="response",s=cv$lambda.min)#predict
  plot(fit)
  library(plotmo) # for plotres
  plotres(fit)
  dev.off()
  #par(mar=c(4.5,4.5,1,4))
  plot(fit)
  vnat=coef(fit)
  vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
  vn=paste("var",1:24)
  axis(4, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5)
  pred
target
  result.difference<-pred-target
  result.difference
  #result.difference  


cv
ggplot
lattice , coefficients, 
predicted-target -gene




plot(pred[,1], target)
abline(lm( target ~ pred[,1]), col=2)
summary(lm( target ~ pred[,1]))
cor.test(pred[,1], target)
pred[1,]
target-pred

opt.lam = c(cv$lambda.min, cv$lambda.1se)
coef(cv, s = opt.lam)


