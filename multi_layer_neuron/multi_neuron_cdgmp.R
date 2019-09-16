snp.data<-read.csv("cdg-data/cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
#snp.data3<-read.csv("cdg-data/cdgSNPmatrix-Jinyuan.csv3", sep =",", header = T, row.names = 1)
t.snp.data<-t(snp.data)
cdgmp.data<-read.csv("cdg-data/cdgTable.csv2", sep =",", header =T)
gene1 <- read.csv("cdg-data/gene-1.snp.csv", header = T, row.names = 1)

boxplot(cdgmp.data[,2] ~ cdgmp.data[,1], las =2)
abline(h=c(-1.5,-.5), col=2)

m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])

m.cdgmp<-m.cdgmp[row.names(t.snp.data)]
length(m.cdgmp[row.names(t.snp.data)])
rownames(t.snp.data)

# h_dgmp<-m.cdgmp[m.cdgmp > -.5]
# m_dgmp<-m.cdgmp[m.cdgmp < -.5 & m.cdgmp > -1.5]
# l_dgmp<-m.cdgmp[m.cdgmp < -1.5 ]

t1 <- as.matrix(ifelse(m.cdgmp > -.5, 1,0 ))
t2 <- as.matrix(ifelse(m.cdgmp < -.5 & m.cdgmp > -1.5, 1,0))
t3 <- as.matrix(ifelse(m.cdgmp  < -1.5 , 1,0 ))

gene1 <- read.csv("cdg-data/gene-1.snp.csv2", header = T, row.names = 1)
snps <- t(gene1)


rownames(snps)
colnames(snps)
ncol(snps)

w <- runif(32, 1e-3, 1e-2)
b <- runif(1)
weight_bias<-as.matrix(c(w,b))

snps2<-cbind(snps,rep(1,30))

linear.combination<-snps2 %*% weight_bias
linear.combination





rm(list=ls(all=TRUE))
x<-iris
eta=0.1#the learning rate
targets=data.frame(t1=c(rep(1,50),rep(0,50),rep(0,50)),t2=c(rep(0,50),rep(1,50),rep(0,50)),t3=c(rep(0,50),rep(0,50),rep(1,50)))
gneuron<-function(n){
  w<-runif(4, 1e-3, 1e-2)#random weights
  b<-runif(1)#random bias
  x1<-iris[,1];x2<-iris[,2]; x3<-iris[,3]; x4<-iris[,4]
  z=x1*w[1] + x2*w[2] + x3*w[3] + x4*w[4] +b
  h<-list(z,w,b)
  return(h)
}

l.neurons<-lapply(1:3, function(x) gneuron(x));names(l.neurons)<-c("z1","z2","z3")

weight_bias<-lapply(l.neurons, function(dx){
  c(unlist(dx[2]),unlist(dx[3]))
})

weight_bias
m.weight_bias=matrix(unlist(weight_bias), ncol =3)
m.weight_bias
x2<-as.matrix(x[,1:4], ncol =4);x2<-cbind(x2,rep(1, 150))

m.weight_bias=matrix(unlist(weight_bias), ncol =3)
wdf=data.frame(w1=numeric(), w2=numeric(),w3=numeric(),w4=numeric(), b=numeric()) 
nwdf<-list(wdf,wdf,wdf)
accuracy<-list()
x2
m.weight_bias
y
targets
y
#for (generation in 1:1000){
  linear.combination<-x2 %*% m.weight_bias  
  y <- exp(linear.combination) / rowSums(exp(linear.combination))
  y
  e<-targets-y
  for (neuron in 1:3){
    m.weight_bias[5,neuron]=m.weight_bias[5,neuron]-eta*-sum(e[,neuron])/150    
    for (weight in 1:4){
      gc<-x2[,weight]*e[,neuron]
      m.weight_bias[weight,neuron]=m.weight_bias[weight,neuron] - eta* (-sum(gc)/150)
    }
    nwdf[[neuron]]<-rbind(nwdf[[neuron]],t(m.weight_bias)[neuron,])    
  }  
  accuracy[[generation]]<-length(which(max.col(y)==max.col(targets)))/150
}


# y
# as.matrix(nwdf[[1]][1:5][1000,])
# x2 %*% m.weight_bias
# targets
# 
# 
# sp0<-x[which(x$Species == "setosa"),]
# sp1<-x[which(x$Species == "versicolor"),]
# sp2<-x[which(x$Species == "virginica"),]
# sp0


















# for (i in length(m.cdgmp)){
#   if (m.cdgmp > -.5){
#   array[i]<-1
#   }
#   else if (m.cdgmp < -.5 & m.cdgmp > -1.5){
#   }
# }
# t <- as.matrix(ifelse(cdg.mean[row.names(snps)] < -1, 0,1 ))
# 
# df<-data.frame(m_dgmp)
# 
# mf<-as.matrix(h_dgmp, ncol =1)
# 
# 
# t.snp.data
# dim(t.snp.data)
# 
# 
# colnames(t.snp.data)
# 
# ?gsub
# 
# gsub("?<|fig||", "", colnames(t.snp.data),perl=T)#, fixed = FALSE, perl = FALSE)
# t.snp.data
# #<-m.cdgmp[m.cdgmp[,2] > -.5,1:2]
# #m_dgmp<-m.cdgmp[m.cdgmp < -.5 & m.cdgmp > -1.5,1:2]
# 
# 
# 
# 
# 
# 
# h_dgmp
# snp.gene1<-read.csv("cdg-data/gene-1.snp.csv2", sep =",", header = T, row.names = 1)
# t.snp.gene1<-t(snp.gene1)
# t.snp.gene1
# m.cdgmp
# t.snp.data
# #gene-1.snp.csv2
# #t <- as.matrix(ifelse(cdg.mean[row.names(snps)] < -1, 0,1 ))
