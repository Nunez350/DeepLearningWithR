rm(list=ls(all=TRUE))

library(glmnet)
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp_data.R')
tg<-t(snp.by.gene.uniq[[49]])#2 dim 3 by 30
t<-m.cdgmp
x<-tg#for testing purposes
c.lasso<-function(x){
exit
  eix
split<-sample(nrow(x), floor(0.7*nrow(x)))
train.d<-as.matrix(cbind(x,t)[split,])
test.d<-as.matrix(cbind(x,t)[-split,])
train.m<-train.d[,-ncol(train.d)]
target<-train.d[,ncol(train.d)]
target<-abs(target)
fit <- glmnet(train.m, target, family = "poisson")
cv<-cv.glmnet(train.m,target,nfolds = 10)#cross_validate to get the best lambda
pred<-predict(fit, train.m, type="response",s=cv$lambda.min)#predict
return(list(fit=fit,cv=cv,pred=pred ))
}
plot(cv)
output<-lapply(snp.by.gene.uniq, function(x){
  ac<-c.lasso(t(x))
  print(ac)
})
str(output)

par(mfrow=c(2,5))  
par(mar=c(.1,.1,.1,.1))
par(mai=c(8,8))
par(maa)
output$fit
lasso_plots<-function(x,s){
  lapply(snp.by.gene.uniq, function(x){
    x<-t(x)
    

    i=parent.frame()$i[]
    split<-sample(nrow(x), floor(0.7*nrow(x)))
    train.d<-as.matrix(cbind(x,t)[split,])
    test.d<-as.matrix(cbind(x,t)[-split,])
        train.m<-train.d[,-ncol(train.d)]
    target<-train.d[,ncol(train.d)]
    target<-abs(target)
    fit <- glmnet(train.m, target, family = "poisson")
    cv<-cv.glmnet(train.m,target,nfolds = 10)#cross_validate to get the best lambda
    pred<-predict(fit, train.m, type="response",s=cv$lambda.min)#p
    
#    snps<-nrow(x)[i]
    #    snps<-nrow(x[i])
    print(ncol(x))
    plot(fit, label = TRUE) 
    title(paste("Gene",i,"SNPS",ncol(x)),line = +2)

    #return(plot(fit, label = TRUE, main=c(paste("gene",i),paste("SNPS",ncol(x)))))
#    return(plot(fit, label = TRUE, main=c(paste("Gene",i,"SNPS",ncol(x) ))))
        #return(plot(fit, label = TRUE) , title( main=c(paste("Gene",i,"SNPS",ncol(x) )), line = +1))
    #title( main=c(paste("Gene",i,"SNPS",ncol(x) )), line = +1)
#    return(plot(fit, label = TRUE, main=c(paste("gene",i),paste("SNPS",snps)))
    i<-i+1
})

}


par(mfrow=c(2,5))  
#par(mai=c(1,1))  
lasso_plots(snp.by.gene.uniq)    
dev.off()

return(plot(
      
      tsne(d.gene,initial_dims=nrow(x), perplexity = 2), main=c(paste("gene ",i), 
                                                                          paste("SNPS",nrow(x)))

lapply()
cv
fit
plot(fit$beta,ylab=expression(beta),)
plot(fit, label = TRUE)

plot(fit$a0)
#summary(fit)

fit$beta

plot(cv)
plot(fit, xvar = "lambda", label = TRUE)
plot(fit, xvar = "dev", label = TRUE)



library(plotmo) # for plotres
plotres(fit)
dev.off()
#par(mar=c(4.5,4.5,1,4))

plot(pred[,1], target)
abline(lm( target ~ pred[,1]), col=2)
summary(lm( target ~ pred[,1]))
r<-cor.test(pred[,1], target)
r$parameter
r$p.value
r$estimate


pred[,1]
pred

vnat
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
# 
# 
# lm_corr_plot<-function(x,s){
#   lapply(snp.by.gene.uniq, function(x){
#     x<-t(x)
#     
#     i=parent.frame()$i[]
#     split<-sample(nrow(x), floor(0.7*nrow(x)))
#     train.d<-as.matrix(cbind(x,t)[split,])
#     test.d<-as.matrix(cbind(x,t)[-split,])
#     train.m<-train.d[,-ncol(train.d)]
#     target<-train.d[,ncol(train.d)]
#     target<-abs(target)
#     fit <- glmnet(train.m, target, family = "poisson")
#     cv<-cv.glmnet(train.m,target,nfolds = 10)#cross_validate to get the best lambda
#     pred<-predict(fit, train.m, type="response",s=cv$lambda.min)
#     plot(pred[,1], target)
#     tryCatch(
#       title(paste("Gene",i,"SNPS",ncol(x)),line = +2)
#       abline(lm( target ~ pred[,1]), col=2)
#       ,err=function(e) NULL
#     )
#     #summary(lm( target ~ pred[,1]))
#     #r<-cor.test(pred[,1], target)
#     #r$parameter
#     #r$p.value
#     #r$estimate
#     
#     #plot(fit, label = TRUE) 
#     #title(paste("Gene",i,"SNPS",ncol(x)),line = +2)
#     i<-i+1
#     
#   })
#   
# }
# lm_corr_plot(snp.by.gene.uniq)