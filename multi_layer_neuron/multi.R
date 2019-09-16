rm(list=ls(all=TRUE))
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
#library(gsubfn)
# #get_snps_targets<-function(){
#snps<-list()
#targets<-list()
get_snps_targets<-function(){
snp.data<-read.csv("cdg-data/cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
t.snp.data<-t(snp.data)
cdgmp.data<-read.csv("cdg-data/cdgTable.csv2", sep =",", header =T)
gene1 <- read.csv("cdg-data/gene-1.snp.csv", header = T, row.names = 1)

m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[row.names(t.snp.data)]

t1 <- as.matrix(ifelse(m.cdgmp > -.5, 1,0 ))
t2 <- as.matrix(ifelse(m.cdgmp < -.5 & m.cdgmp > -1.5, 1,0))
t3 <- as.matrix(ifelse(m.cdgmp  < -1.5 , 1,0 ))


gene1 <- read.csv("cdg-data/gene-1.snp.csv2", header = T, row.names = 1)
snps <- t(gene1)
targets<-cbind(t1,t2,t3)
return(list(snps, targets))
}

list[snps,targets]<-get_snps_targets()
snps
targets
#with(get_snps_targets(), a )
#attach(get_snps_targets())


w <- runif(32, 1e-3, 1e-2)
b <- runif(1)
weight_bias<-as.matrix(c(w,b))

snps2<-cbind(snps,rep(1,30))

linear.combination<-snps2 %*% weight_bias
