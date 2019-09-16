rm(list=ls(all=TRUE))
setwd("~/msk-snp-cdgmp/cdg-data/")

snp.data<-read.csv("cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
t.snp.data<-t(snp.data)
snp.dist<-dist(t.snp.data, method= "manhattan")
dim(unique(snp.data[1:32,]))
dim(snp.data[1:32,])

catch<-NULL


  i = 1
  unique(snp.data[i:31,])
dim(snp.data)
rownames(snp.data)
regexpr("fig|287.2645.peg.1474|C576T",rownames(snp.data))
?grep
catch
#gene1 <- read.csv("snp-mei.csv", header = F, row.names = 1)
unq.snps <- read.csv("snp-mei.csv", header = F, row.names = 1)
header <- read.csv("header.csv", header = F)
colnames(unq.snps) <- as.character(as.matrix(header[2:31]))
snps<-t(unq.snps)

# dd<-dist(snps, method = 'manhattan')
# hdd<-hclust(dd)
# plot(hdd)
# snp.dist<-dist(snps, method= "manhattan")
# 
# 
# plot(snp.dist)
# df<-hclust(snp.dist)
# 
# plot(df)
# ct<-cutree(df, k=4)
# length(which(ct ==4))
# length(snps[which(ct == 1)])

cdgmp.data<-read.csv("cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
#excluding cdgmp strains that do not have snp sites
m.cdgmp<-m.cdgmp[row.names(snps)]
sample(rep(1:2, c(20,10))
#randomly splitting the SNPS strains into training and predictive testing sets
split.indeces<-split(1:30, sample(rep(1:2, c(20,10))))
split.indeces
training.snps<-snps[split.indeces$`1`,] # N =20
test.snps<-snps[split.indeces$`2`,]
rownames(snps)
rownames(training.snps)
training.cdgmp<-m.cdgmp[row.names(training.snps)]
testing.cdgmp<-m.cdgmp[row.names(test.snps)]

#reducing the data dimensionality, grouping by hclust
dist.training.snps<-dist(training.snps, method = 'manhattan')
cluster.dist.training.snps<-hclust(dist.training.snps)
plot(cluster.dist.training.snps)
ct<-cutree(cluster.dist.training.snps, k=4)
length(which(ct ==1))
length(snps[which(ct == 1)])



#declaring targets for di-gmp levels high, medium low
t1 <- as.matrix(ifelse(training.cdgmp > -.5, 1,0 )) # high
t2 <- as.matrix(ifelse(training.cdgmp < -.5 & training.cdgmp > -1.5, 1,0)) # medium 
t3 <- as.matrix(ifelse(training.cdgmp  < -1.5 , 1,0 )) #low
targets<-cbind(t1,t2,t3)

eta=0.1#the learning rate
weight_bias<-NULL
for(i in 1:3){
  w <- runif(1078, 1e-3, 1e-2)
  b <- runif(1)
  weight_bias<-cbind(weight_bias,as.matrix(c(w,b)))
}
training.snps2<-cbind(training.snps,rep(1,20))


wdf=data.frame(w1=numeric(), w2=numeric(),w3=numeric(),w4=numeric(), b=numeric()) 
nwdf<-list(wdf,wdf,wdf)
accuracy<-list()
#training.snps2 %*% weight_bias
#dim(training.snps2 %*% weight_bias)
#dim(weight_bias)
#dim(training.snps2)
for (generation in 1:1000){
  linear.combination<-training.snps2 %*% weight_bias
  
  y <- exp(linear.combination) / rowSums(exp(linear.combination))
  e<-targets-y
  for (neuron in 1:3){
    weight_bias[5,neuron]=weight_bias[5,neuron]-eta*(-sum(e[,neuron])/20)
    for (weight in 1:1078){
      gc<-training.snps2[,weight]*e[,neuron]
      weight_bias[weight,neuron]=weight_bias[weight,neuron] - eta* (-sum(gc)/20)
    }
    nwdf[[neuron]]<-rbind(nwdf[[neuron]],t(weight_bias)[neuron,])    
  }  
  accuracy[[generation]]<-length(which(max.col(y)==max.col(targets)))/20
}

plot(0:0, xlim = c(1,1000), ylim=c(-.1.5,.1.5) ,ylab= "weights", xlab="Epochs", main = "Weight Updates")
lapply(nwdf[[1]][1:4],function(x){ lines(1:1000, x, col="2" ) })
lapply(nwdf[[2]][1:4],function(x){ lines(1:1000, x, col="3" ) })
lapply(nwdf[[3]][1:4],function(x){ lines(1:1000, x, col="4") })
legend("topleft", legend = c("neuron 1","neuron 2", "neuron 3"), col=c("2","3","4"), pch=10, horiz = T, cex = .9) 

plot(1:1000, unlist(accuracy), xlab = "Epochs", ylab = "Accuracy", main = "Accuracy Plot")

test.snps<-cbind(test.snps,rep(1,10))

dim(test.snps)
dim(weight_bias)
linear.combination.predict<-test.snps %*% weight_bias
testing.cdgmp
y <- exp(linear.combination.predict) / rowSums(exp(linear.combination.predict))


test1 <- as.matrix(ifelse(testing.cdgmp > -.5, 1,0 )) # high
test2 <- as.matrix(ifelse(testing.cdgmp < -.5 & testing.cdgmp > -1.5, 1,0)) # medium 
test3 <- as.matrix(ifelse(testing.cdgmp  < -1.5 , 1,0 )) #low
test.targets<-cbind(test1,test2,test3)


round(y)
test.targets
length(which(max.col(y)==max.col(test.targets)))/20
dim(snps)

tt<-rbind(snps[1,1:10] ,snps[2,1:6],snps[5,1:6])
tt<-snps[1:5,1:10]

unique(t(tt))
dim(unique(t(snps)))

dim(snps)
which(tt unique()
save<-NULL
for ( i in ncol(tt)){
  save[i]<-unique(tt[,i])  
  
}
save

save<-list()

unique(snps[,1])
unique(snps[,2])

Grp <- vapply(unique(df$G), function(x) paste(df$SKU[which(df$G==x)], collapse = ""), "abc", USE.NAMES = FALSE)
vapply(unique())
?vapply

out<-apply(snps[,1:1078],2, function(xx){
  unique(xx)
})
out[1,1:10]
?apply
dim(out)
for ( i in 1:4){
  save[i]<-unique(snps[,i])  
  
}
save
weight_bias[1000,]

#using weight and bias to predict species 
l.bias<- as.numeric(bias.list[1000])
fw1<-as.numeric(weights1[1000]); fw2<-as.numeric(weights2[1000])
y=1/(1+exp(-(fw1*x1 + fw2 * x2 + l.bias)))


round(y,1)
##   [1] 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 1 1
##  [36] 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
##  [71] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1


sample(rep(1:2, c(20,10)))


#scratch code 1
#snp.data<-read.csv("cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
snp.data<-read.csv("cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)


x="fig|287.2645.peg.2872|C642T"

nn<-rownames(snp.data)
nn
u.names<-unique(strtrim(nn, 8))
u.names
class(unique(strtrim(nn, 8)))



f<-factor(1:51)
length(u.names)
#substring(x,14,21)
?factor
f<-factor(x= u.names)
f
nn
for ( i in 1:10){
  
  #output<-
  output<-sapply(snp.data, function(xx){
    
    output2<-unique(strtrim(rownames(xx), 8))#rownames(snp.data)))
    f<-factor(output2)
    return(f)
  })
  f
  sapply( f, function(xx){
    snp.data[f]
    
  })
  output
  output2
  unique(strtrim(rownames(snp.data),8))
  
  
}

?hash
library(hash)
strtrim(rownames(snp.data)[1], 8)
unique(strtrim(rownames(snp.data), 8))




j

collect<-list()
c2<-list()
gene<-strtrim(rownames(snp.data)[1], 8) 
for ( j in 1:3) {
  gene=temp
  for ( i in 1:100) {
    count<-0
    if (strtrim(rownames(snp.data)[i],8)== gene   )  {
      print(c("i-",i,"j-",j ))
      print(strtrim(rownames(snp.data)[i],8))
      collect[count]<-strtrim(rownames(snp.data)[i],8)
      count<-count+1
    } 
    print("no")
    else{
      temp =strtrim(rownames(snp.data)[i],8)
      
    }
    
  }
}

collect
strtrim(rownames(snp.data), 8)


table(collect)

strtrim(rownames(snp.data), 8)
+}








