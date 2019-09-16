rm(list=ls(all=TRUE))
setwd("~/Users/roynunez/cleanNow/msk-snp-cdgmp/cdg-data/")
snp.data<-read.csv("/Users/roynunez/cleanNow/msk-snp-cdgmp/cdg-data/cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
test1 <- snp.by.gene.uniq[[1]]
test1<-t(test1)

cdgmp.data<-read.csv("/Users/roynunez/cleanNow/msk-snp-cdgmp/cdg-data/cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(test1)]


eta=0.1#the learning rate
#based on genes
snp.by.gene.uniq
##dist(snp.by.gene.uniq)

# #reducing the data dimensionality, grouping by hclust
#dist.training.snps<-dist(training.snps, method = 'manhattan')
# cluster.dist.training.snps<-hclust(dist.training.snps)
# plot(cluster.dist.training.snps)
# ct<-cutree(cluster.dist.training.snps, k=4)
# length(which(ct ==1))
# length(snps[which(ct == 1)])

test1 <- snp.by.gene.uniq[[1]]
snps<-t(test1)

neuralFunction<-function(x){
  performance<-list()
  hit.acc<-list()
  test.hit.acc<-list()
  for (trials in 1:10) {
    
    #randomly splitting the SNPS strains into training and predictive testing sets
    split.indeces<-split(1:30, sample(rep(1:2, c(20,10))))
    training.snps<-snps[split.indeces$`1`,] # N =20
    test.snps<-snps[split.indeces$`2`,]
    training.cdgmp<-m.cdgmp[rownames(training.snps)]
    testing.cdgmp<-m.cdgmp[rownames(test.snps)]
    
    #declaring targets for di-gmp levels high, medium low
    t1 <- as.matrix(ifelse(training.cdgmp > -.5, 1,0 )) # high
    t2 <- as.matrix(ifelse(training.cdgmp < -.5 & training.cdgmp > -1.5, 1,0)) # medium 
    t3 <- as.matrix(ifelse(training.cdgmp  < -1.5 , 1,0 )) #low
    targets<-cbind(t1,t2,t3)
    
    weight_bias<-NULL
    for(i in 1:3){
      w <- runif(ncol(training.snps), 1e-3, 1e-2)
      b <- runif(1)
      weight_bias<-cbind(weight_bias,as.matrix(c(w,b)))
    }
    
    
    training.snps2<-cbind(training.snps,rep(1,20))
    
    wdf=data.frame(w1=numeric(), w2=numeric(),w3=numeric(),w4=numeric(), b=numeric()) 
    nwdf<-list(wdf,wdf,wdf)
    accuracy<-list()
    training.snps2  
    weight_bias
    for (generation in 1:1000){
      linear.combination<-training.snps2 %*% weight_bias
      
      y <- exp(linear.combination) / rowSums(exp(linear.combination))
      e<-targets-y
      for (neuron in 1:3){
        weight_bias[5,neuron]=weight_bias[5,neuron]-eta*(-sum(e[,neuron])/20)
        for (weight in 1:ncol(training.snps)){
          gc<-training.snps2[,weight]*e[,neuron]
          weight_bias[weight,neuron]=weight_bias[weight,neuron] - eta* (-sum(gc)/20)
        }
        nwdf[[neuron]]<-rbind(nwdf[[neuron]],t(weight_bias)[neuron,])    
      }  
      accuracy[[generation]]<-length(which(max.col(y)==max.col(targets)))/20
    }
    
    
    
    final.wts <- matrix(c(nwdf[[1]][1000,], nwdf[[2]][1000,], nwdf[[3]][1000,]), nrow=7, byrow=F)
    
    final.wts2 <- matrix(unlist(final.wts), nrow=7, byrow = F)
    L<-training.snps2 %*% final.wts2
    y <- exp(L) / rowSums(exp(L))
    hits<-length(which(round(y[,1])==targets[,1] & round(y[,2])==targets[,2] & round(y[,3])==targets[,3]))
    hit.acc<-rbind(hit.acc, hits/20)
    
    
    
    test1 <- as.matrix(ifelse(testing.cdgmp > -.5, 1,0 )) # high
    test2 <- as.matrix(ifelse(testing.cdgmp < -.5 & testing.cdgmp > -1.5, 1,0)) # medium 
    test3 <- as.matrix(ifelse(testing.cdgmp  < -1.5 , 1,0 )) #low
    test.targets<-cbind(test1,test2,test3)
    
    test.snps<-cbind(test.snps,rep(1,10))
    
    test.line<-test.snps %*% final.wts2
    test.y <- exp(test.line) / rowSums(exp(test.line))
    test.hits<-length(which(round(test.y[,1])==test.snps[,1] & round(test.y[,2])==test.snps[,2] & round(test.y[,3])==test.snps[,3]))
    
    test.hit.acc<-rbind(test.hit.acc,test.hits/10, deparse.level = 1)
  }
  performance<-list(hit.acc, range(hit.acc), test.hit.acc, range(test.hit.acc))
  names(performance)<-c("self-prediction","range_self_predict", "external.prediction", "range_external_prediction")
  return(performance)
}


out<-lapply(snp.by.gene.uniq, function(x)
  neuralFunction(snps)
)

out.df <- lapply(out, function(x) data.frame(self=unlist(x$`self-prediction`[,1]), target=unlist(x$`external.prediction`[,1])))

for(i in 1:length(out.df)) {
  n <- names(out.df[i]);
  boxplot(out.df[[i]], ylim=c(0,1), ylab="accuracy", las=1, main=n);
  stripchart(out.df[[i]], vertical = T, method = 'jitter', add = TRUE, pch=16, col=1:2, ylim=c(0,1));
}






library(ape)
setwd("~/msk-snp-cdgmp/cdg-data/")
tr=read.tree("cdg-tree-v1-mid.dnd")
txt.names=read.table("cdg.strains.txt3", sep="\t", header = F, row.names = 1)
tr$tip.label <- as.character(txt.names[tr$tip.label,]) #
#tr$tip.label<-as.character(txt.names$V2[match(tr$tip.label, txt.names$V1)])
plot(tr, font =1)
add.scale.bar()

library(phylobase)
library(phylotools)
tr$tip.label<-gsub("pa_.+_", "",tr$tip.label)
t.snp.data<-t(snp.data)
ace(t.snp.data[tr$tip.label,1],tr,type= "d")

plot(tr)
column=1
co <- c("blue", "yellow")
tr.unrooted<-unroot(tr)
plot.mpr <- function(column=1) {
  plot.phylo(tr, main = colnames(t.snp.data)[column])
  tmpr<-MPR(t.snp.data[,column], tr.unrooted, outgroup = "F34365")
  #nodelabels(paste("[", tmpr[, 1], ",", tmpr[, 2], "]", sep = "")) 
  tiplabels(t.snp.data[,column][tr.unrooted$tip.label], adj = -2)
  return(tmpr)
}





# 
# tmpr
# class(tr$tip.label)
# dim(tr$tip.label)
# match(tr$tip.label, t.snp.data)
# snp.data
# 
# rownames(t.snp.data)
# tr$tip.label
# match(rownames(t.snp.data), tr$tip.label)
# length(match(rownames(t.snp.data), tr$tip.label))
# ace([tr.unrooted$tip.label,1],tr,type= "d")
# snp.data[tr$tip.label,1]
# dim(t.snp.data)
# 
# tr$tip.label
# 
# class(t.snp.data[,1])
# dim(t.snp.data[,1])
# structure(t.snp.data[,1])
# 
# tr$tip.label
# dim(t(snp.data))
# ######ace(t.snp.data[tr$tip.label,1],tr,type= "d")
# 
# dim(as.matrix(t(snp.by.gene.uniq)))
# ace(t(snp.data,tr, type="d"))
# tr
# dim(t(snp.data))
# ace([tr.unrooted$tip.label,1],tr,type= "d")
# length(snp.by.gene.uniq)
# tr
# as.matrix(tr)
# ph.t <- t(ph.mat)
# rownames(ph.t) <- gsub("X", "x", rownames(ph.t))
# rownames(t.snp.data)[tr$tip.label,1]
# dim(t.snp.data[,1])
# 
# ace(t.snp.data[tr$tip.label,1],tr,type="d")
# t.snp.data[tr$tip.label,1]
# strep.ace<-ace(ph.t[tr.unrooted$tip.label,1],tr,type= "d")
# co <- c("blue", "yellow")
