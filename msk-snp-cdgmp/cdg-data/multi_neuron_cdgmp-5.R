rm(list=ls(all=TRUE))
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
setwd("~/msk-snp-cdgmp/cdg-data/")

get_snps_targets<-function(){
  snp.data<-read.csv("cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
  t.snp.data<-t(snp.data)
  cdgmp.data<-read.csv("cdgTable.csv2", sep =",", header =T)
  gene1 <- read.csv("gene-1.snp.csv", header = T, row.names = 1)
  
  m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
  names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
  m.cdgmp<-m.cdgmp[row.names(t.snp.data)]
  
  t1 <- as.matrix(ifelse(m.cdgmp > -.5, 1,0 ))
  t2 <- as.matrix(ifelse(m.cdgmp < -.5 & m.cdgmp > -1.5, 1,0))
  t3 <- as.matrix(ifelse(m.cdgmp  < -1.5 , 1,0 ))
  
  
  gene1 <- read.csv("gene-1.snp.csv2", header = T, row.names = 1)
  snps <- t(gene1)
  targets<-cbind(t1,t2,t3)
  return(list(snps, targets))
}

list[snps,targets]<-get_snps_targets()

weight_bias<-NULL
for(i in 1:3){
  w <- runif(32, 1e-3, 1e-2)
  b <- runif(1)
  weight_bias<-cbind(weight_bias,as.matrix(c(w,b)))
}

snps2<-cbind(snps,rep(1,30))

eta=0.1#the learning rate
wdf=data.frame(w1=numeric(), w2=numeric(),w3=numeric(),w4=numeric(), b=numeric()) 
nwdf<-list(wdf,wdf,wdf)
accuracy<-list()
for (generation in 1:1000){
  linear.combination<-snps2 %*% weight_bias
  y <- exp(linear.combination) / rowSums(exp(linear.combination))
  e<-targets-y
  for (neuron in 1:3){
    weight_bias[5,neuron]=weight_bias[5,neuron]-eta*(-sum(e[,neuron])/30)
    for (weight in 1:32){
      gc<-snps2[,weight]*e[,neuron]
      weight_bias[weight,neuron]=weight_bias[weight,neuron] - eta* (-sum(gc)/30)
    }
    nwdf[[neuron]]<-rbind(nwdf[[neuron]],t(weight_bias)[neuron,])    
  }  
  accuracy[[generation]]<-length(which(max.col(y)==max.col(targets)))/30
}

plot(0:0, xlim = c(1,1000), ylim=c(-3,3.5) ,ylab= "weights", xlab="Epochs", main = "Weight Updates")
lapply(nwdf[[1]][1:4],function(x){ lines(1:1000, x, col="2" ) })
lapply(nwdf[[2]][1:4],function(x){ lines(1:1000, x, col="3" ) })
lapply(nwdf[[3]][1:4],function(x){ lines(1:1000, x, col="4") })
legend("topleft", legend = c("neuron 1","neuron 2", "neuron 3"), col=c("2","3","4"), pch=10, horiz = T, cex = .9) 

plot(1:1000, unlist(accuracy), xlab = "Epochs", ylab = "Accuracy", main = "Accuracy Plot")

