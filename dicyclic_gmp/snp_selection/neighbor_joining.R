rm(list=ls(all=TRUE))
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp.R')
library(ape)
gene1<-t(snp.by.gene.uniq[[1]])
gene1<-t(gene1)
gene1
tr=nj(dist.gene(gene1))
gene1
tr
plot(tr)
min(tr$edge.length)
tr$edge.length
which(tr$edge.length == min(tr$edge.length))
m.outgroup<-tr$tip.label[which(tr$edge.length == min(tr$edge.length))]
m.outgroup
r.tr<-root(tr,m.outgroup,resolve.root = TRUE)
r.tr
#o.tr<-order(tr$edge.length)
o.tr
#na.omit(tr$tip.label[o.tr])
f=multi2di(r.tr, random = TRUE)
plot(f)
scale()
par(mfrow=c(1,1))
library(phangorn)
?midpoint()


replicate(10, rTraitDisc(tr, states = c(0,1), rate = 100, model = "ER"))
