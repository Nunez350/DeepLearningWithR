rm(list=ls(all=TRUE))

library("tsne")
#library("Rtsne")
snp.data<-read.csv("~/msk-snp-cdgmp/cdg-data/cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
snps<-lapply(snp.by.gene.uniq, function(x) t(x))
cdgmp.data<-read.csv("~/msk-snp-cdgmp/cdg-data/cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(snps[[1]])]
gene1<-snp.by.gene.uniq[[1]]
groups <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))

mt<-t(gene1)

split<-sample(nrow(mt), floor(0.70*nrow(mt)))

train<-mt[split,]
test<-mt[-split,]
dim(train)
#?tsne
d.gene1<-dist(snp.by.gene.uniq[[1]])
gene50<-snp.by.gene.uniq[[50]]
gene1<-snp.by.gene.uniq[[1]]
d.gene50<-dist(snp.by.gene.uniq[[50]])
n1<-31
n2<-40



par(mfrow=c(2,5))  

tsne_plots<-function(x,s,e){
  lapply(x[s:e], function(x){
    i=parent.frame()$i[]
    d.gene<-dist(t(x))
    i<-i+(s-1)
    return(plot(tsne(d.gene,initial_dims=nrow(x), perplexity = 2), main=c(paste("gene ",i), 
                                                                          paste("SNPS",nrow(x)))
                
    ))
    #return(plot(tsne(d.gene)))
  }
  )
}  

run.tsne.plot<-function(){
  pdf("TSNE_SNPS")
  par(mfrow=c(2,5)) 
  tsne_plots(snp.by.gene.uniq, 1,10)
  tsne_plots(snp.by.gene.uniq, 11,21)
  tsne_plots(snp.by.gene.uniq, 22,33)
  tsne_plots(snp.by.gene.uniq, 34,45)
  #par(mfrow=c(1,5))
  tsne_plots(snp.by.gene.uniq, 46,50)
  dev.off()
  par(mfrow=c(1,1))
}
#run.tsne.plot()
