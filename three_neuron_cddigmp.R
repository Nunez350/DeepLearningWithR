rm(list=ls(all=TRUE))
setwd("~/msk-snp-cdgmp/cdg-data/")
snp.data<-read.csv("cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
test1 <- snp.by.gene.uniq[[1]]
test1<-t(test1)

cdgmp.data<-read.csv("cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(test1)]


eta=0.1#the learning rate




# #reducing the data dimensionality, grouping by hclust
# dist.training.snps<-dist(training.snps, method = 'manhattan')
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
    split.indeces
    snps
    training.snps<-snps[split.indeces$`1`,] # N =20
    training.snps
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
}

lapply(out.df,function(dx){
  n <- names(dx);
  #boxplot(dx, ylim=c(0,1), ylab="accuracy", las=1, main=n);
#stripchart(dx[[1]], vertical = T, method = 'jitter', add = TRUE, pch=16, col=1:2, ylim=c(0,1))
})

#pdf("Self prediction vs External Prediction")
lapply(out,function(dx){
  plot(dx[[1]], type="o", col="blue", ylim=c(0,1), axes=FALSE, ann=FALSE)
  lines(dx[[3]], type="o", pch=22, lty=2, col="red")
  title(main="predictions", col.main="red", font.main=4)
  title(xlab="percent accuracy", col.lab=rgb(0,0.5,0))
  title(ylab="Iterations", col.lab=rgb(0,0.5,0))
  legend(1, g_range[2], c("self-prediction","external.prediction"), cex=0.8, 
         col=c("blue","red"), pch=21:22, lty=1:2);
  
})

lapply(out,function(dx){
  g_range<-range(0, dx[[1]],dx[[3]])
  plot(dx[[1]], type="o", col="blue", ylim=g_range, axes=FALSE, ann=FALSE)
  axis(2, las=1, at=4*0:g_range[2])
  lines(dx[[3]], type="o", pch=22, lty=2, col="red")
  title(main="predictions", col.main="red", font.main=4)
  title(xlab="percent accuracy", col.lab=rgb(0,0.5,0))
  title(ylab="Iterations", col.lab=rgb(0,0.5,0))
  legend(1, g_range[2], c("self-prediction","external.prediction"), cex=0.8, 
         col=c("blue","red"), pch=21:22, lty=1:2);
  
})
#dev.off()


library(ape)
tr=read.tree("cdg-tree-v1.dnd")
txt.names=read.table("cdg.strains.txt3")
tr$tip.label<-as.character(txt.names$V2[match(tips, txt.names$V1)])
plot(tr)
# 
# #return(list(snps, targets))
# gene1<-neuralFunction(snps)
# gene1
# plot(gene)
# gene1$self-prediction
# plot(gene1$`self-prediction`, gene1$external.prediction, type= "l")
# lapply(nwdf[[1]][1:4],function(x){ lines(1:1000, x, col="2" ) })
# 
# 
# 
# g_range<-range(0, gene1$`self-prediction`,gene1$external.prediction)
# 
# plot(gene1$`self-prediction`, type="o", col="blue", ylim=g_range, 
#      axes=FALSE, ann=FALSE)
# axis(2, las=1, at=4*0:g_range[2])
# box()
# lines(gene1$external.prediction, type="o", pch=22, lty=2, col="red")
# title(main="predictions", col.main="red", font.main=4)
# title(xlab="percent accuracy", col.lab=rgb(0,0.5,0))
# title(ylab="Iterations", col.lab=rgb(0,0.5,0))
# 
# 
# legend(1, g_range[2], c("self-prediction","external.prediction"), cex=0.8, 
#        col=c("blue","red"), pch=21:22, lty=1:2);
# 
# 
# 
# 
# out
# out$peg.1075
# length(out)
# lappfor (i in 1:50 ){
#   print(names(out[i]))
#   
# }
lapply(out, function(dx){
   names(dx)
 })
# 
# names(out)
# lapply(out,function(dx){
#   g_range<-range(0, dx[[1]],dx[[3]])
#   plot(dx[[1]], type="o", col="blue", ylim=g_range, 
#        axes=FALSE, ann=FALSE)
#   axis(2, las=1, at=4*0:g_range[2])
#   lines(dx[[3]], type="o", pch=22, lty=2, col="red")
#   title(main="predictions", col.main="red", font.main=4)
#   title(xlab="percent accuracy", col.lab=rgb(0,0.5,0))
#   title(ylab="Iterations", col.lab=rgb(0,0.5,0))
#   legend(1, g_range[2], c("self-prediction","external.prediction"), cex=0.8, 
#          col=c("blue","red"), pch=21:22, lty=1:2);
# 
# })
# 
# 
# 
# 
# 
# g_range<-range(0, gene1$`self-prediction`,gene1$external.prediction)
# 
# plot(gene1$`self-prediction`, type="o", col="blue", ylim=g_range, 
#      axes=FALSE, ann=FALSE)
# axis(2, las=1, at=4*0:g_range[2])
# box()
# lines(gene1$external.prediction, type="o", pch=22, lty=2, col="red")
# title(main="predictions", col.main="red", font.main=4)
# title(xlab="percent accuracy", col.lab=rgb(0,0.5,0))
# title(ylab="Iterations", col.lab=rgb(0,0.5,0))
# 
# 
# legend(1, g_range[2], c("self-prediction","external.prediction"), cex=0.8, 
#        col=c("blue","red"), pch=21:22, lty=1:2);
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# test3 <- as.matrix(ifelse(testing.cdgmp  < -1.5 , 1,0 )) 
# check<-matrix()
# check<-list()
# for ( j in 1:10){
# check<-rbind(check,snps[j,1], deparse.level = 1)
# }
# test.hit.acc<-rbind(test.hit.acc,test.hits/10, deparse.level = 1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# linear.combination.predict<-training.snps2 %*% weight_bias
# y <- exp(linear.combination.predict) / rowSums(exp(linear.combination.predict))
# 
# round(y)
# 
# which(round(y[,1])==targets[,1] & round(y[,2])==targets[,2] & round(y[,3])==targets[,3])
# 
# est.tawhich(round(y[,1])==test.targets[,1] & round(y[,2])==test.targets[,2] & round(y[,3])==test.targets[,3] )
# plot(0:0, xlim = c(1,1000), ylim=c(-4,4) ,ylab= "weights", xlab="Epochs", main = "Weight Updates")
# lapply(nwdf[[1]][1:4],function(x){ lines(1:1000, x, col="2" ) })
# lapply(nwdf[[2]][1:4],function(x){ lines(1:1000, x, col="3" ) })
# lapply(nwdf[[3]][1:4],function(x){ lines(1:1000, x, col="4") })
# legend("topleft", legend = c("neuron 1","neuron 2", "neuron 3"), col=c("2","3","4"), pch=10, horiz = T, cex = .9) 
# 
# plot(1:1000, unlist(accuracy), xlab = "Epochs", ylab = "Accuracy", main = "Accuracy Plot")
# 
# test1 <- as.matrix(ifelse(testing.cdgmp > -.5, 1,0 )) # high
# test2 <- as.matrix(ifelse(testing.cdgmp < -.5 & testing.cdgmp > -1.5, 1,0)) # medium 
# test3 <- as.matrix(ifelse(testing.cdgmp  < -1.5 , 1,0 )) #low
# test.targets<-cbind(test1,test2,test3)
# 
# test.snps<-cbind(test.snps,rep(1,10))
# 
# linear.combination.predict<-test.snps %*% weight_biaslinear.combination.predict<-test.snps %*% weight_bias
# y <- exp(linear.combination.predict) / rowSums(exp(linear.combination.predict))
# targets
# 
# 
# 
# 
# #length(which(max.col(y)==max.col(test.targets)))/10
# which(round(y)==(test.targets))
# length(which(round(y)==(test.targets)))/length(test.targets)
# 
# test.targets
# round(y)
# max.col(y)
# max.col(test.targets)
# # 
# # t.snp.data<-t(snp.data)
# # snp.dist<-dist(t.snp.data, method= "manhattan")
# # 
# # 
# # 
# # snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
# # 
# # test <- snp.by.gene.uniq[[1]]
# # test
# 
# 
# final.wts <- matrix(c(nwdf[[1]][1000,], nwdf[[2]][1000,], nwdf[[3]][1000,])), nrow=7, byrow=F)
# 
# final.wts2 <- matrix(unlist(final.wts), nrow=7, byrow = F)
# L<-training.snps2 %*% final.wts2
# y <- exp(L) / rowSums(exp(L))
# which(round(y[,1])==targets[,1] & round(y[,2])==targets[,2] & round(y[,3])==targets[,3])
