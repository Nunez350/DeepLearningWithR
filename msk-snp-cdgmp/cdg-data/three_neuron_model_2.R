rm(list=ls(all=TRUE))
setwd("~/msk-snp-cdgmp/cdg-data/")
snp.data<-read.csv("cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
snps<-lapply(snp.by.gene.uniq, function(x) t(x))
head(snps)
snps
cdgmp.data<-read.csv("cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(snps[[1]])]
m.cdgmp
gene1<-snp.by.gene.uniq[[1]]
gene50<-snp.by.gene.uniq[[50]]
eta=0.1#the learning rate

weight_bias_list<-NULL
weight_bias_list<- lapply(snps, function(x){
  split.indeces<-split(1:30, sample(rep(1:2, c(20,10))))
  training.snps<-x[split.indeces$`1`,]
  #training.snps<-x[[1]][split.indeces$`1`,]
  training.snps<-cbind(training.snps,co=rep(1,ncol(training.snps)))
  training.snps
  #training.snps
  test.snps<-x[split.indeces$`2`,]
  #test.snps<-x[[1]][split.indeces$`2`,]
  test.snps<-cbind(test.snps,co=rep(1,nrow(test.snps)))
  training.cdgmp<-m.cdgmp[rownames(training.snps)]
  testing.cdgmp<-m.cdgmp[rownames(test.snps)]  
  t1 <- as.matrix(ifelse(training.cdgmp > -.5, 1,0 )) # high
  t2 <- as.matrix(ifelse(training.cdgmp < -.5 & training.cdgmp > -1.5, 1,0)) # medium 
  t3 <- as.matrix(ifelse(training.cdgmp  < -1.5 , 1,0 )) #low
  targets<-cbind(t1,t2,t3)
  training.snps
  neuron1 <- runif(ncol(training.snps)-1, 1e-3, 1e-2)
  neuron2 <- runif(ncol(training.snps)-1, 1e-3, 1e-2)
  neuron3 <- runif(ncol(training.snps)-1, 1e-3, 1e-2)
  neuron1<- append(neuron1,runif(1))
  neuron2<- append(neuron2,runif(1))
  neuron3<- append(neuron3,runif(1))
  weight_bias<-cbind(neuron1,neuron2, neuron3)
  #  training.snps  
  
  return(list(targets=targets,training.snps=training.snps, crossvvalidatesnps=test.snps,weight_bias=weight_bias))
})
#gene=1;weight=1; neuron=1
#
gene=1
as.matrix(weight_bias_list[[gene]]$training.snps)%*% as.matrix(weight_bias_list[[gene]]$weight_bias)
ff=    weight_bias_list[[1]]$weight_bias
length(weight_bias_list[1:2])
for (gene in 1:length(weight_bias_list[1:2])){
  for (generation in 1:10000){
    nweights=ncol(weight_bias_list[[gene]]$training.snps)-1
    linear.combination<- as.matrix(weight_bias_list[[gene]]$training.snps)%*% as.matrix(weight_bias_list[[gene]]$weight_bias)
    y <- exp(linear.combination) / rowSums(exp(linear.combination))
    e<-weight_bias_list[[1]]$targets-y
    
    for (neuron in 1:3){
      weight_bias_list[[gene]]$weight_bias[nweights+1,neuron]=weight_bias_list[[gene]]$weight_bias[nweights+1,neuron]-eta*(-sum(e[,neuron])/20)
      
      
      for (weight in 1:nweights){
        gc<-  weight_bias_list[[gene]]$weight_bias[weight,neuron]*e[,neuron]
        weight_bias_list[[gene]]$weight_bias[weight,neuron]=weight_bias_list[[gene]]$weight_bias[weight,neuron]-eta*(-sum(gc)/20)
      }
    }
  }  
}
ff-    weight_bias_list[[1]]$weight_bias

weight_bias_list[[1]]
weight_bias_list[[1]]$weight_bias
weight_bias_list[[1]]$crossvvalidatesnps
weight_bias_list[[1]]$training.snps
dim(weight_bias_list[[3]]$training.snps)
gene
stor<-list()
for (gene in 1:50){
L<-weight_bias_list[[gene]]$training.snps %*%weight_bias_list[[gene]]$weight_bias
#weight_bias_list[[1]]$crossvvalidatesnps%*%weight_bias_list[[1]]$weight_bias
#weight_bias_list[[1]]$targets
y <- exp(L) / rowSums(exp(L))

rowSums(exp(L))
stor[[gene]]=length(which(round(y[,1])==weight_bias_list[[gene]]$targets[,1] & round(y[,2])==weight_bias_list[[gene]]$targets[,2] & round(y[,3])==weight_bias_list[[gene]]$targets[,3]))/nrow(weight_bias_list[[gene]]$targets)

}

stor
plot()
stor
?max.col
length(which(round(y[,1])==weight_bias_list[[1]]$targets[,1] & round(y[,2])==weight_bias_list[[1]]$targets[,2] & round(y[,3])==weight_bias_list[[1]]$targets[,3]))/nrow(weight_bias_list[[1]]$targets)

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
  neuralFunction(x)
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
# 
# 
# ##dist(snp.by.gene.uniq)

# #reducing the data dimensionality, grouping by hclust
#dist.training.snps<-dist(training.snps, method = 'manhattan')
# cluster.dist.training.snps<-hclust(dist.training.snps)
# plot(cluster.dist.training.snps)
# ct<-cutree(cluster.dist.training.snps, k=4)
# length(which(ct ==1))
# length(snps[which(ct == 1)])

#neuralFunction<-function(x){
performance<-list()
# dim(weight_bias_list[[1]]$training.snps)
# weight_bias_list[[1]]$training.snps
# dim(weight_bias_list[[1]]$weight_bias)
# weight_bias_list[[1]]$weight_bias
# 
# as.matrix(weight_bias_list[[1]]$training.snps)%*% as.matrix(weight_bias_list[[1]]$weight_bias)

# dim(as.matrix(weight_bias_list[[1]]$training.snps))
# as.matrix(weight_bias_list[[1]]$training.snps)
# dim(as.matrix(weight_bias_list[[1]]$weight_bias))
# as.matrix(weight_bias_list[[1]]$weight_bias)
