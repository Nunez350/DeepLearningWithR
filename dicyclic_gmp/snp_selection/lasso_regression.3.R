#lasso
rm(list=ls(all=TRUE))
library(glmnet)
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp.R')
t <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))
class(t)
dim(t)

t<-as.matrix(m.cdgmp)
tg<-t(snp.by.gene.uniq[[49]])
#dim(snp.by.gene.uniq[[49]])
x<-tg

c.lasso<-function(x){
  x
  dim(x)
  hist(t(t))
  plot(t)
  dim(t)
  split<-sample(nrow(x), floor(0.7*nrow(x)))
  train.d<-as.matrix(cbind(x,t)[split,])
  dim(train.d)
  test.d<-as.matrix(cbind(x,t)[-split,])
  dim(train.d)
  train.d[,ncol(train.d)]
  ?glmnet
  fit2 <- glmnet(train.d, train.d[,ncol(train.d)], family = "gaussian")
  fit2 <- glmnet(train.d, train.d[,ncol(train.d)], family = "poisson")
  dim(train.d)
  summary(fit2)
  cv.2<-cv.glmnet(train.d,train.d[,ncol(train.d)],nfolds = 10)#cross_validate to get the best lambda
  plot(cv.2)
  summary(cv.2)
  
  pred<-predict(fit2, train.d, type="response",s=cv.2$lambda.min)#predict
  pred
  
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
fit2$betaplot(1:90, fit2$beta[1,], type = "l", ylim=c(-9,8))
fit2$beta
lout
fit2$a0






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


plot(1:90, fit2$beta[1,], type = "l", ylim=c(-9,8))
lines(1:90, fit2$beta[2,], type = "l", col=2)
lines(1:90, fit2$beta[3,], type = "l", col=3)
lines(1:90, fit2$beta[4,], type = "l", col=4)
lines(1:90, fit2$beta[5,], type = "l", col=5)
lines(1:90, fit2$beta[6,], type = "l", col=6)
