p.value(gene1, df =1)

library("tsne")
library("Rtsne")
?tsne
d.gene1<-dist(snp.by.gene.uniq[[1]])
gene50<-snp.by.gene.uniq[[50]]
gene1<-snp.by.gene.uniq[[1]]
d.gene50<-dist(snp.by.gene.uniq[[50]])
par(mfrow=c(1,1))
n1<-31
n2<-40

test<-lapply(snp.by.gene.uniq[n1:n2], function(x){

    d.gene<-dist(t(x))
    return(plot(tsne(d.gene,initial_dims=nrow(x), perplexity = 2), main=n1))
  #return(plot(tsne(d.gene)))
}
)

test<-lapply(snp.by.gene.uniq[n1:n2], function(x){
  
  d.gene<-dist(t(x))
  return(plot(princomp(t(gene1))$scores))
  #return(plot(tsne(d.gene)))
}
)
gene1

plot(test,main="1-10")
par(mfrow=c(2,5))  


plot(tsne(t(gene1)))
_______________
snp.by.gene.uniq[1:10]
lapply(snp.by.gene.uniq, function(x){
  #dim(snp.by.gene.uniq)
  dim(x)
})


?tsne

d.gene<-dist(gene1)
plot(tsne(d.gene,initial_dims=50, perplexity = 30))
length((d.gene[[1]]))
nrow(d)
dim(gene1)
gene1

## Not run: 
colors = rainbow(length(unique(iris$Species)))
names(colors) = unique(iris$Species)
ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }

tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)

# compare to PCA
dev.new()
#pca_iris = 
princomp(iris[,1:4])  #$scores[,1:2]
princomp(t(gene1))$scores[,1:2]
plot(princomp(t(gene1))$scores)

gene1
plot(pca_iris, t='n')
  
text(pca_iris, labels=iris$Species,col=colors[iris$Species])
dev.off()
## End(Not run)
## 
## 
## 
## 
library("Rtsne")
?Rtsne

Rtsne(t(gene1))
dim(gene1)
u.gene1<-unique(gene1[ , ] )
dim(u.gene1)
dim(gene1)
d.gene1<-dist(snp.by.gene.uniq[[1]])
gene50<-snp.by.gene.uniq[[50]]
gene1<-snp.by.gene.uniq[[1]]
d.gene50<-dist(snp.by.gene.uniq[[50]])
dim(unique(gene1[,1:ncol(gene1)]))

dim(snp.by.gene.uniq[[1]])
Rtsne(d.gene1,N=6, is_distance =T,check_duplicates = FALSE,verbose=T, pca_scale=T)
Rtsne(gene1,N=6, is_distance =F,check_duplicates = FALSE, verbose=T, pca_scale=T)

Rtsne(X, dims = 2, initial_dims = 50, perplexity = 30,
      theta = 0.5, check_duplicates = TRUE, pca = TRUE, max_iter = 1000,
      verbose = FALSE, is_distance = FALSE, Y_init = NULL,
      pca_center = TRUE, pca_scale = FALSE,
      stop_lying_iter = ifelse(is.null(Y_init), 250L, 0L),
      mom_switch_iter = ifelse(is.null(Y_init), 250L, 0L), momentum = 0.5,
      final_momentum = 0.8, eta = 200, exaggeration_factor = 12, ...)


iris_unique <- unique(iris) # Remove duplicates
dim(iris_unique)
dim(iris_unique)
iris_matrix <- as.matrix(iris_unique[,1:4])
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(iris_matrix) # Run TSNE
head(tsne_out)
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=iris_unique$Species)

# Using a dist object
tsne_out <- Rtsne(dist(iris_matrix))
plot(tsne_out$Y,col=iris_unique$Species)

# Use a given initialization of the locations of the points
tsne_part1 <- Rtsne(iris_unique[,1:4], theta=0.0, pca=FALSE,max_iter=350)
tsne_part2 <- Rtsne(iris_unique[,1:4], theta=0.0, pca=FALSE, max_iter=150,Y_init=tsne_part1$Y)
m.gene1<-as.matrix(gene1)
Rtsne(m.gene1[,1:5])
m.gene1[,1:10]
?dist
d.gene1<-dist(gene1, )
Rtsne(d.gene1, check_duplicates = F, is_distance =T) 



---------------------
install.packages(c("devtools","broom","MASS"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","snpStats","DESeq2"))
biocLite(c("apply", "anyMissing", "rowMedians"))
library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)
library("SnpMatrix")
data(for.exercise)
m<-signature(from = "t(gene1)", to = "SnpMatrix")
m
class(t(gene1))
tg<-t(gene1)
dim(tg)
y<-as.raw(tg)
snp.gene1<-new("SnpMatrix",y, rownames(y))
?as.raw
y=as.raw(t(gene1))
ibsDist(tg)
class(as.array(tg))
?as.array
ibs.stats(gene1)
View(snp.gene1)
snp.gene1[1,]
xxt(snp.gene1, correct.for.missing = F)
?xxt
eigen(t(gene1), symmetric=TRUE)
?xxt
gene1
use <- seq(1, ncol(snps.10), 10)

length(use)
str(snps.10)


sub.10 <- snps.10[,use]
xxmat <- xxt(t(gene1), correct.for.missing=FALSE)
View(xxmat)
str(sub.10)

xxmat <- xxt(sub.10, correct.for.missing=FALSE)
xxmat
evv <- eigen(xxmat, symmetric=TRUE)
evv
pcs <- evv$vectors[,1:5]
pcs  
library(ape)
t=nj(dist.gene(gene1))##############
t$
par(m)  
snpdata = sub.10@.Data
snpdata
subject.support$cc
status = subject.support$cc
status
snp1 = as.numeric(snpdata[,1])
snp1
snp1[snp1==0] = NA
dim(m.cdgmp)
dim(gene1)
glm1 = glm(status ~ snp1,family="binomial")
glm1 = glm(t ~ tg,family="binomial")
glm1
plot(glm1)

tidy(glm1)
m.cdgmp
t <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))
tg
t
tf<-factor(t)

library("glmnet")
tf<-as.factor(t)
df<-cbind(df,t)
df<-as.data.frame(cbind(tg,t))
df
sparse.model.matrix(~.,df, sparsem=T)
library(Matrix)
mt<-cbind(tg,t)

set.seed(2)
split<-sample(nrow(mt), floor(0.7*nrow(mt)))
split
train<-mt[split,]
test<-mt[-split,]
class(train)
train<-df[split,]
test<-df[-split,]
dim(train)
dim(test)
sparse.model.matrix(~.,train, sparsem=T)
sparse.model.matrix(~.,test[1:10])
library(glmnet)
as.matrix(train)
mt<-as.matrix(train)
dim(mt)
dim(mtest)
head(mtest)
mtest<-as.matrix(test)

fit<-glmnet(mt, mt[,7])######
fit
cv<-cv.glmnet(mt,mt[,7],nfolds = 3)

pred<-predict(fit, mtest, type="response",s=cv$lambda.min)
round(pred)==mtest[,7]
library(pROC)
plot(mtest[,7],pred)
roc(mtest[,7],pred)
class(as.matrix(mtest[,7]))

class(train[,7])
train
glmnet()
rownames(df)
glmnet()
cv<-cv.glmnet(tg, t, alpha = 1, nlambda = 100)#mean squared error
cv<-cv.glmnet(tg, t, alpha = 1, family="binomial",type.measure= "class",nlambda = 100)
plot(cv)
fit<-glmnet(tg, t, alpha = 1, lambda = cv$lambda.1se)
fit<-glmnet(tg, t, alpha = 1, family="binomial",lambda =  cv$lambda.1se)

fit$beta[,1]
plot(fit$beta)
`plot(fit)
cv$lambda.1se


ridge.mod<-glmnet(tg, tf, alpha = 1, nlambda = 100)
tf
ridge.mod
predict(ridge.mod, s = 0, exact = T, type = 'coefficients')

cdgmp.data$logcdg
status
snp1
tidy(glm1)
}
[
  }
)
#Random Forest
can be used for classification or regression
avoids overfitting
can deal with a large number of features 
helps with feature selection based on importance
#mtry predictors at each node
#
# default sq.root(p)
#p/3 for regression
#
data<-c(1,2,3,4,5)
#
#NSP is the response variable
as.factor(data$NSP)
table(data$NSP)
ind <_sampe(2,nrow(data), repalce = T,prob = c(0.7, 0.3))
train<-data[ind==1,]
test<-data[ind==2,]
#1!ntree bottstrap samples
#for eacn bootstrap sample 