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
m.weight_bias=matrix(unlist(weight_bias), ncol =3)

x2<-as.matrix(x[,1:4], ncol =4);x2<-cbind(x2,rep(1, 150))

m.weight_bias=matrix(unlist(weight_bias), ncol =3)
wdf=data.frame(w1=numeric(), w2=numeric(),w3=numeric(),w4=numeric(), b=numeric()) 
nwdf<-list(wdf,wdf,wdf)
accuracy<-list()
for (generation in 1:1000){
  linear.combination<-x2 %*% m.weight_bias  
  y <- exp(linear.combination) / rowSums(exp(linear.combination))
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
iris


plot(0:0, xlim = c(1,1000), ylim=c(-3,3.5) ,ylab= "weights", xlab="Epochs")
lapply(nwdf[[1]][1:4],function(x){ lines(1:1000, x, col="2", lty = sample(1:4, replace = F)) })
lapply(nwdf[[2]][1:4],function(x){ lines(1:1000, x, col="3", lty = sample(1:4, replace = F)) })
lapply(nwdf[[3]][1:4],function(x){ lines(1:1000, x, col="4", lty = sample(1:4, replace = F)) })
legend("topleft", legend = c("neuron 1","neuron 2", "neuron 3"), col=c("2","3","4"), pch=10, horiz = T, cex = .9) 

plot(1:1000, unlist(accuracy))

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
















