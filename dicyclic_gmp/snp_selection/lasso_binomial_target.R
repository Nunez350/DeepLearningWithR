#lasso for binary target
rm(list=ls(all=TRUE))
library(glmnet)
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp_data.R')
t <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))
tg<-t(snp.by.gene.uniq[[49]])#2 dim 3 by 30
x<-tg#for testing purposes
x[,25]
dim(x)
c.lasso<-function(x){
  split<-sample(nrow(x), floor(0.7*nrow(x)))
  train.d<-as.matrix(cbind(x,t)[split,])
  test.d<-as.matrix(cbind(x,t)[-split,])
  fit <- glmnet(train.d, train.d[,ncol(train.d)], family = "binomial")
  plot(fit$a0)
  fit$npasses
  plot(fit$beta)
  plot(fit)
  plot(train.d[,-ncol(train.d)])
  fit$beta[1,]
  fit$beta
  plot(1:100, fit$beta[1,], type = "l", ylim=c(-9,8))
  length(fit$beta[1,])
  cv<-cv.glmnet(train.m,target,nfolds = 10)#cross_validate to get the best lambda
  plot(cv.2)
  plot(cv)
  
  # After running cv.glmnet, you don't have to rerun glmnet! Every lambda in the grid (cv$lambda) has already been run. This technique is called "Warm Start" and you can read more about it here. Paraphrasing from the introduction, the Warm Start technique reduces running time of iterative methods by using the solution of a different optimization problem (e.g., glmnet with a larger lambda) as the starting value for a later optimization problem (e.g., glmnet with a smaller lambda).
  summary(cv)
  plot(fit)  
  plot(fit, xvar = "lambda", label = TRUE)
  # Let’s plot “fit” against the log-lambda value and with each curve labeled.
  plot(fit, xvar = "dev", label = TRUE)
  # Now when we plot against %deviance we get a very different picture. This is percent deviance explained on
  # the training data. What we see here is that toward the end of the path this value are not changing much, but
  # the coefficients are “blowing up” a bit. This lets us focus attention on the parts of the fit that matter. This
  # will especially be true for other models, such as logistic regression.
  # 
  
  #cv$lambda.1se instead of cv$lambda.min
  pred<-predict(fit, train.d, type="response",s=cv.2$lambda.min)#predict
  pred
  predict.fishnet(fit, train.d, type="response",s=cv.2$lambda.min)#predict
  print(fit)
  plot(fit)
  train.d[,ncol(train.d)]
  pred==train.d[,ncol(train.d)]
  round(pred)==train.d[,ncol(train.d)]
  mod.prediction<-round(pred)==train.d[,ncol(train.d)]
  length(which(mod.prediction==TRUE))
  
  acc<-length(round(pred)==nrow(train.d))/nrow(train.d)*100
return(acc)
}
c.lasso(x)
dim(fit2$beta)
plot(1:length(o),o$acc1)
plot(1:100, fit2$beta[1,], type = "l", ylim=c(-9,8))
lout<-list()
par(mfrow=c(1,1)) 
summary(fit2)
lapply(snp.by.gene.uniq, function(x){
  ac<-c.lasso(t(x))
  print(ac)
})
fit2$beta
plot(1:90, fit2$beta[1,], type = "l", ylim=c(-9,8))
fit2$beta
lout
fit$a0



pred<-predict(fit, train.m, type="response",s=cv.2$lambda.min)#predict
pred

dim(train.m)
library(SimPhe)

x1 <- rnorm(4000, mean = 5, sd = 10)
x2 <- rnorm(4000, mean = 10, sd = 30)
x <- matrix(cbind(x1, x2), ncol = 2)
# test original correlation
cor.test(x[, 1], x[, 2])
# correlation matrix
corM <- matrix(c(1, 0.6, 0.6, 1), ncol = 2)
# standard deviation matrix
sdM <- matrix(c(10, 0, 0, 30), ncol = 2)
# build correlation
x.new <- build.cor.phe(x, corM, sdM)
# check mean and standard deviation of new data set
apply(x.new, 2, mean)
apply(x.new, 2, sd)
# test correlation
cor.test(x.new[, 1], x.new[, 2])
x.new[, 1]

x1 <- rnorm(4000, mean = 5, sd = 10)
x2 <- rnorm(4000, mean = 10, sd = 30)
x <- matrix(cbind(x1, x2), ncol = 2)
build.sd.matrix(x)


plot(xx, yy, type = "l")  ## plot density curve (or use `plot(d)`)

plot
plot(fit)
plot(glmnet(x,y),label = TRUE )
plot(1:90, fit2$beta[1,], type = "l", ylim=c(-9,8))
lines(1:90, fit2$beta[2,], type = "l", col=2)
lines(1:90, fit2$beta[3,], type = "l", col=3)
lines(1:90, fit2$beta[4,], type = "l", col=4)
lines(1:90, fit2$beta[5,], type = "l", col=5)
lines(1:90, fit2$beta[6,], type = "l", col=6)
data.matrix(l=fit$lambda, d=deviance(fit))[abs(l-reg.lambda) == min(abs(l-reg.lambda))]$d[1]
l=fit$lambda
l
#d=deviance(fit))[abs(l-reg.lambda) 
df<-matrix(l,d)
d=deviance(fit)
df[abs(l-reg.lambda) == min(abs(l-reg.lambda))]$d[1]

library(plotmo) # for plotres
plotres(fit)
par(mar=c(4.5,4.5,1,4))
plot(fit)
vnat=coef(fit)
vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
vn=paste("var",1:24)
axis(4, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5)
length(vnat)
Options are almost the same as the Gaussian family except that for type.measure, * “deviance” (default)
gives the deviance * “mse” stands for mean squared error * “mae” is for mean absolut

opt.lam = c(cv.2$lambda.min, cv.2$lambda.1se)
coef(cv.2, s = opt.lam)
glmnet.control()
#lmnet fits the entire solution path for Lasso or elastic-net problems efficiently with various
techniques such as warm starts, which are essential for the nonlinear models like logistic regression