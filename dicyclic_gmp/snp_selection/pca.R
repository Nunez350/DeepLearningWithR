rm(list=ls(all=TRUE))
source('~/machine_learning/dicyclic_gmp/snp_selection/di_ci_gmp.R')

test_pca_plots<-function(x,s,e){
  #  par(mfrow=c(2,5)) 
  lapply(snp.by.gene.uniq[s:e], function(x){
    i=parent.frame()$i[]
    i<-i+(s-1)
    #d.gene<-dist(t(x))
    #if(princomp(t(x))$scores){
    pca_scores<-princomp(t(x))$scores
    #}else{te<-prcomp(t(x))$scores}
    
    #pca_scores<-princomp(t(x))$scores
    
    return(tryCatch(plot(pca_scores,main=c(paste("gene",i),paste("SNPS",nrow(x)))), error=function(e) NULL))    
  }
  )
}

?princomp()  

test_pca_plots(snp.by.gene.uniq,1,9)

#return(plot(tsne(d.gene)))
s<-1;e<-3
pca_scores<-list()
test_pca_plots<-function(x,s,e){
  #  par(mfrow=c(2,5)) 
  out<-lapply(x, function(x){
    i=parent.frame()$i[]
    i<-i+(s-1)
    if(nrow(x) > ncol(x)){
      princomp(x)$scores
    }
    else{return(("NA"))}
    #return(tryCatch(plot(pca_scores,main=c(paste("gene",i),paste("SNPS",nrow(x)))), error=function(e) NULL))    
  }
  )
  return(out)
}
nn<-test_pca_plots(t(snp.by.gene.uniq),1,50)

nn[1:9]
t(snp.by.gene.uniq[[1]])

pca_scores<-test_pca_plots(snp.by.gene.uniq,1,50)
if()
  tryCatch(tt!=0)
dim(snp.by.gene.uniq[[8]])
try(princomp(t(snp.by.gene.uniq[[8]])))
princomp(t(snp.by.gene.uniq[[9]]))
?try
?tryCatch()

dim(pca_scores)
str(pca_scores)

par(mfrow=c(2,5)) 
test_pca_plots<-function(x,s,e){

  lapply(x[s:e], function(x){
  i=parent.frame()$i[]
  i<-i+(s-1)
  d.gene<-dist(t(x))
  return(plot(princomp(t(x))$scores,
              main=c(paste("gene ",i), paste("SNPS",nrow(x))
                     )
              ))
  #return(plot(tsne(d.gene)))
    }
  }
)
par(mfrow=c(2,5)) 
n1
n2
x<-gene1
t(x)
test_pca_plots<-function(x,s,e){
#  par(mfrow=c(2,5)) 
  lapply(snp.by.gene.uniq[s:e], function(x){
  i=parent.frame()$i[]
  i<-i+(s-1)
  #d.gene<-dist(t(x))
  if(princomp(t(x))$scores){
    pca_scores<-princomp(t(x))$scores
  }else{te<-prcomp(t(x))$scores}
  
  #pca_scores<-princomp(t(x))$scores
  
  
  return(plot(pca_scores,main=c(paste("gene",i),paste("SNPS",nrow(x))))
  }#else{next}
    #return(plot(tsne(d.gene)))
)
}

?princomp()  
test_pca_plots(snp.by.gene.uniq,1,10)
  par(mfrow=c(1,1)) 

  if(princomp(t(x))$scores){
    te<-princomp(t(x))$scores
  }else{te<-prcomp(t(x))$scores}

  plot(prcomp(t(snp.by.gene.uniq))$scores)
dim(te)
dim(princomp(t(x))$scores)
prcomp(t(x))$scores

princomp(t(snp.by.gene.uniq[[1]]))$scores
prcomp(snp.by.gene.uniq[[9]])$scores

if(exists(princomp(t(snp.by.gene.uniq[[1]]))$scores)){
  print(prprincomp(t(snp.by.gene.uniq[[1]]))$scores)
}


plot(princomp(snp.by.gene.uniq[[8]])$scores)
plot(princomp(snp.by.gene.uniq[[9]])$scores)
dim(snp.by.gene.uniq[[9]])
dim(gene1)
test
tt<-princomp(snp.by.gene.uniq[[9]])$scores
if(exists("tt")){
  print(tt)
}else{next}

plot(test,main="1-10")
plot(1:1000,1:1000)
par(mfrow=c(1,1))  






# test<-function(x,s,e){
#   lapply(x[s:e], function(x){
#     i=parent.frame()$i[]
#     d.gene<-dist(t(x))
#     i<-i+(s-1)
#     return(plot(tsne(d.gene,initial_dims=nrow(x), perplexity = 2), main=c(paste("gene ",i), 
#                                                                           paste("SNPS",nrow(x)))
#                 
#     ))
#     #return(plot(tsne(d.gene)))
#   }
#   )
# } 
# test(snp.by.gene.uniq,1,10)

exit()
<_
<-
  