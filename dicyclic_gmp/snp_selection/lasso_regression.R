#lasso
rm(list=ls(all=TRUE))
library(glmnet)
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp.R')
t <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))

dim(snp.by.gene.uniq[[49]])
tg<-t(snp.by.gene.uniq[[49]])
x<-tg

c.lasso<-function(x){
  x
  split<-sample(nrow(x), floor(0.7*nrow(x)))
  train.d<-as.matrix(cbind(x,t)[split,])
  test.d<-as.matrix(cbind(x,t)[-split,])
  fit2 <- glmnet(train.d, train.d[,ncol(train.d)], family = "binomial")
  cv.2<-cv.glmnet(train.d,as.matrix(train.d[,ncol(train.d)]),nfolds = 10)#cross_validate to get the best lambda
  pred<-predict(fit2, train.d, type="response",s=cv.2$lambda.min)#predict
  round(pred)==train.d[,ncol(train.d)]
  cv.2<-cv.glmnet(train.d[,-ncol(train.d)],train.d[,ncol(train.d)],nfolds = 10)#cross_validate to get the best lambda
#fit.2<-glmnet(train.d[,-ncol(train.d)], train.d[,ncol(train.d)], family = "binomial")


#fit<-glmnet(mtrain, mtrain[,ncol(mtrain)])######
#cv<-cv.glmnet(mtrain,mtrain[,ncol(mtrain)],nfolds = 10)#cross_validate to get the best lambda
#dim(mtrain)

d#im(mtrain[,-ncol(mtrain)])

plot(cv)
pred<-predict(fit, mtest, type="response",s=cv$lambda.min)#predict
pred<-predict(fit2, train.d, type="response",s=cv.2$lambda.min)#predict
fit2
dim(train.d)
pred
plot(cv)
pred
fit
summary(fit)
round(pred)==mtest[,-ncol(mtrain)]
mtest[,ncol(mtrain)]
prediction<-round(pred)==mtest[,ncol(mtrain)]
accuracy<-length(which(prediction == TRUE))/length(prediction)*100
pred2<-predict(fit, mtrain, type="response",s=cv$lambda.min)#predict
round(pred2)==mtrain[,-ncol(mtrain)]
prediction2<-round(pred2)==mtrain[,ncol(mtrain)]
accuracy2<-length(which(prediction2 == TRUE))/length(prediction2)*100


t.m<-as.matrix(t(gene1))
t.m%*%t
co<-cor(t.m[,1:6])

which(co > .6, arr.ind = T)

dim(t.m)
dim(t)
return(list(acc1=accuracy, acc2=accuracy2))
}
c.lasso(t(snp.by.gene.uniq[[1]]))
t(snp.by.gene.uniq[[1]])
getwd()

snp.by.gene.uniq
h<-list()
for (i in 1:1000){
split.indeces<-sample(length(snp.by.gene.uniq), floor(0.7*length(snp.by.gene.uniq)))
trial1<-snp.by.gene.uniq[split.indeces]
trial1.2<-snp.by.gene.uniq[-split.indeces]

out<-lapply(trial1, function(gene_snp){
    x<-t(gene_snp)
    c.lasso(x)
})
a1<-lapply(out,function(dx){
  return(dx$acc1)
})
h[[i]]<-a1

}
h
table(unlist(h))

lapply(out, function(o){
  plot(1:length(o),o$acc1)
})


x<-NULL
o<-out
plot(out)
a1<-lapply(out,function(dx){
  return(dx$acc1)
})

a2<-lapply(out,function(dx){
  return(dx$acc2)
})
plot(1:35,unlist(a1))
class(a1)

# #predict against itself
# pred2<-predict(fit, mtrain, type="response",s=cv$lambda.min)
# round(pred2)==mtrain[,ncol(mtrain)]
# 
# r<-as.matrix(cbind(tg,t)[split,])
# pred3<-predict(fit, r, type="response",s=cv$lambda.min)
# round(pred3)==r[,ncol(r)]
# x<-NULL
# snp.by.gene.uniq
# lapply(snp.by.gene.uniq, function(gene_snp){
# gene_snp<-snp.by.gene.uniq[[1]]
#   x<-t(gene_snp)
#   x
#   r<- as.matrix(cbind(x,t))
#   r
#   fit
#   pred3<-predict(fit, r, type="response",s=cv$lambda.min)
#   round(pred3)==r[,ncol(r)]
#   }
# )
# 
# dim(mtest)
# 
# c.lasso(tg)
# t(snp.by.gene.uniq[[8]])
# tg<-t(snp.by.gene.uniq[[8]])
# tg
# save<-lapply(snp.by.gene.uniq, function(x){
#   c.lasso(t(x))
# }
# )
# save
# length(save)
# 
# 
# which(save ==TRUE)
# r<-as.matrix(cbind(tg,t)[split,])
# pred3<-predict(fit, r, type="response",s=cv$lambda.min)
# round(pred3)==r[,ncol(r)]
