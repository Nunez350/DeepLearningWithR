snp.data<-read.csv("cdgSNPmatrix-Jinyuan.csv", sep =",", header = T, row.names = 1)
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





