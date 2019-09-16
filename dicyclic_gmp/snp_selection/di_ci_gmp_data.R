#di-ci-gmp data file
#
#rm(list=ls(all=TRUE))
getwd()
#library("tsne")
#library("Rtsne")
snp.data<-read.csv("~/msk-snp-cdgmp/cdg-data/cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
dim(unlist(snp.by.gene))
str(snp.by.gene)
snps<-lapply(snp.by.gene.uniq, function(x) t(x))
cdgmp.data<-read.csv("~/msk-snp-cdgmp/cdg-data/cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(snps[[1]])]
gene1<-snp.by.gene.uniq[[1]]
groups <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))

mt<-t(gene1)
#library(Matrix)

split<-sample(nrow(mt), floor(0.70*nrow(mt)))

train<-mt[split,]
test<-mt[-split,]
dim(train)
#?tsne
d.gene1<-dist(snp.by.gene.uniq[[1]])
gene50<-snp.by.gene.uniq[[50]]
gene1<-snp.by.gene.uniq[[1]]
d.gene50<-dist(snp.by.gene.uniq[[50]])
#p.value(gene1, df =1)
#n1<-31
#n2<-40
