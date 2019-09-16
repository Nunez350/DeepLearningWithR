rm(list=ls(all=TRUE))
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

h_dgmp<-m.cdgmp[m.cdgmp > -.5]
m_dgmp<-m.cdgmp[m.cdgmp < -.5 & m.cdgmp > -1.5]
l_dgmp<-m.cdgmp[m.cdgmp < -1.5 ]

t1 <- as.matrix(ifelse(m.cdgmp > -.5, 1,0 ))
t2 <- as.matrix(ifelse(m.cdgmp < -.5 & m.cdgmp > -1.5, 1,0))
t3 <- as.matrix(ifelse(m.cdgmp  < -1.5 , 1,0 ))

gene1 <- read.csv("cdg-data/gene-1.snp.csv2", header = T, row.names = 1)

snps <- t(gene1)
rownames(snps)
#colnames(snps)
#row.names(gene1)


ncol(gene1)
nrow(gene1)

row.names(gene1)

w <- runif(32, 1e-3, 1e-2)
b <- runif(1)





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
