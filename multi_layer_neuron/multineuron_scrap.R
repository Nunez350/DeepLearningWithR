rm(list=ls(all=TRUE))
x<-iris;eta=0.1#the learning rate
t1<-c(rep(1,50),rep(0,50),rep(0,50));t2<-c(rep(0,50),rep(1,50),rep(0,50));t3<-c(rep(0,50),rep(0,50),rep(1,50))#third target

#generates random weights, initializes weights, random bias and activation
gneuron<-function(n){
  w<-runif(4, 1e-3, 1e-2)#random weights
  #w1=w[1];w2=w[2];w3=w[3];w4=w[4]#initiating weights
  b<-runif(1)#random bias
  x1<-iris[,1];x2<-iris[,2]; x3<-iris[,3]; x4<-iris[,4]
  z=x1*w[1] + x2*w[2] + x3*w[3] + x4*w[4] +b
  h<-list(z,w,b)
  return(h)
}
l.neurons<-lapply(1:3, function(x) gneuron(x));names(l.neurons)<-c("z1","z2","z3")

l.neurons$z1

for (generation in 1:1000){
  #l.neurons<-lapply(1:3, function(x) gneuron(x));names(l.neurons)<-c("z1","z2","z3")
  z1e<-exp(unlist(l.neurons[[1]][1]))
  z2e<-exp(unlist(l.neurons[[2]][[1]]))
  z3e<-exp(unlist(l.neurons[[3]][[1]]))
  z1e
  exp.act.list<-list(z1e,z2e,z3e)
  
  y1<-unlist(z1e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));
  y2<-unlist(z2e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));
  y3<-unlist(z3e)/(unlist(z1e)+unlist(z2e)+unlist(z3e))
  y=cbind(unlist(y1),unlist(y2),unlist(y3))
  e1=t1-y[,1];e2=t2-y[,2];e3=t3-y[,3];e<-cbind(e1,e2,e3)
  updated.wx<-NULL;gc<-NULL;delta.w<-NULL;n.r<-NULL;sum.w<-NULL;bias<-NULL;delta.bias<-NULL
  for (n in 1:3){
    delta.bias[[n]]<-append(bias,(-sum(e[,n])/150))
    
    for (c in 1:4){
      for (r in 1:150){
        gc<-rbind(gc, x[r,c] * e[r,n])
      }
      delta.w[c]<- -sum(gc)  
      gc<-NULL
    }
    n.r[[n]]<-delta.w
    updated.wx<-NULL
    sw<-NULL
  }
  n.r
  delta.w
  delta.bias
  #w_current - eta * (-sum(t - y) * x)
  #weights_biases[1:nrow(weights_biases) - 1, neuron] - eta * (-colSums(x[, 1:ncol(x) - 1] * e[, neuron])) / nrow(t)
  
  for (n in 1:3){
    l.neurons[[n]][[3]]<-cbind(l.neurons[[n]][[3]], l.neurons[[n]][[3]][generation]-eta*delta.bias[[n]])      
    l.neurons[[n]][[2]]<-rbind(l.neurons[[n]][[2]],l.neurons[[n]][[2]][generation]- eta*n.r[[n]])
    
  }
  
}
l.neurons$z1

plot(1:1001,l.neurons[[1]][[2]][,1], ylim=c(-300,300), col="red", ylab= "weights", xlab="Epochs" )
lines(1:1001,l.neurons[[1]][[2]][,2], col="green")
lines(1:1001,l.neurons[[1]][[2]][,3], col="blue")
lines(1:1001,l.neurons[[1]][[2]][,4], col="orange")
lines(1:1001,l.neurons[[2]][[2]][,1], col="green")
lines(1:1001,l.neurons[[2]][[2]][,2], col="green")
lines(1:1001,l.neurons[[2]][[2]][,3], col="green")
lines(1:1001,l.neurons[[2]][[2]][,4], col="green")
lines(1:1001,l.neurons[[3]][[2]][,1], col="green")
lines(1:1001,l.neurons[[3]][[2]][,2], col="green")
lines(1:1001,l.neurons[[3]][[2]][,3], col="green")
lines(1:1001,l.neurons[[3]][[2]][,4], col="green")
l.neurons[[3]][[2]]




rm(list=ls(all=TRUE))
x<-iris
eta=0.1#the learning rate
t1<-c(rep(1,50),rep(0,50),rep(0,50))
t2<-c(rep(0,50),rep(1,50),rep(0,50))
t3<-c(rep(0,50),rep(0,50),rep(1,50))#third target

#generates random weights, initializes weights, random bias and activation
gneuron<-function(n){
  w<-runif(4, 1e-3, 1e-2)#random weights
  #w1=w[1];w2=w[2];w3=w[3];w4=w[4]#initiating weights
  b<-runif(1)#random bias
  x1<-iris[,1];x2<-iris[,2]; x3<-iris[,3]; x4<-iris[,4]
  z=x1*w[1] + x2*w[2] + x3*w[3] + x4*w[4] +b
  h<-list(z,w,b)
  return(h)
}

l.neurons<-lapply(1:3, function(x) gneuron(x));names(l.neurons)<-c("z1","z2","z3")
l.neurons[[1]][[5]]<-l.neurons[[1]][[2]]
l.neurons[[2]][[5]]<-l.neurons[[2]][[2]]
l.neurons[[3]][[5]]<-l.neurons[[3]][[2]]

for (generation in 1:1000){
  
  for ( n in 1:3){
    temp<-NULL
    for (xw in 1:4){
      temp<-cbind(temp,matrix(l.neurons[[n]][[2]][xw] %*% x[,xw], ncol =1 ))
    }
    l.neurons[[n]][[4]]<-temp[,1]+temp[,2]+temp[,3]+temp[,4]+l.neurons[[n]][[3]]
  }
  
  z1e<-exp(unlist(l.neurons[[1]][[4]]))
  z2e<-exp(unlist(l.neurons[[2]][[4]]))
  z3e<-exp(unlist(l.neurons[[3]][[4]]))
  y1<-unlist(z1e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));
  y2<-unlist(z2e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));
  y3<-unlist(z3e)/(unlist(z1e)+unlist(z2e)+unlist(z3e))
  y=cbind(unlist(y1),unlist(y2),unlist(y3))
  
  e1=t1-y[,1];e2=t2-y[,2];e3=t3-y[,3];e<-cbind(e1,e2,e3)
  updated.wx<-NULL;gc<-NULL;delta.w<-NULL;n.r<-NULL;sum.w<-NULL;bias<-NULL;delta.bias<-NULL
  for (n in 1:3){
    delta.bias[[n]]<-append(bias,(-sum(e[,n])/150))
    for (c in 1:4){
      for (r in 1:150){
        gc<-rbind(gc, x[r,c] * e[r,n])
      }
      delta.w[c]<- -sum(gc)  
      gc<-NULL
    }
    n.r[[n]]<-delta.w
    updated.wx<-NULL
    sw<-NULL
  }
  
  
  for (n in 1:3){
    l.neurons[[n]][[3]]= l.neurons[[n]][[3]]-eta*delta.bias[[n]]
    l.neurons[[n]][[2]]=l.neurons[[n]][[2]]- eta*n.r[[n]]
    
  }
  l.neurons[[n]][[5]]<-rbind(l.neurons[[n]][[5]],l.neurons[[n]][[2]])
}

l.neurons[[n]][[5]]
for(i in 1:10){
  for(n in 1:3){
    
    l.neurons[[n]][[5]]<-rbind(l.neurons[[n]][[5]],l.neurons[[n]][[2]]) 
  }
}





l.neurons[[1]][[5]]


plot(1:1001,l.neurons[[1]][[2]][,1], ylim=c(-300,300), col="red", ylab= "weights", xlab="Epochs" )
lines(1:1001,l.neurons[[1]][[2]][,2], col="green")
lines(1:1001,l.neurons[[1]][[2]][,3], col="blue")
lines(1:1001,l.neurons[[1]][[2]][,4], col="orange")
lines(1:1001,l.neurons[[2]][[2]][,1], col="green")
lines(1:1001,l.neurons[[2]][[2]][,2], col="green")
lines(1:1001,l.neurons[[2]][[2]][,3], col="green")
lines(1:1001,l.neurons[[2]][[2]][,4], col="green")
lines(1:1001,l.neurons[[3]][[2]][,1], col="green")
lines(1:1001,l.neurons[[3]][[2]][,2], col="green")
lines(1:1001,l.neurons[[3]][[2]][,3], col="green")
lines(1:1001,l.neurons[[3]][[2]][,4], col="green")
l.neurons[[3]][[2]]


# l.neurons[[1]][[2]] %*% x[,1:4]
# 
# 
# for ( n in 1:3){
#   temp<-NULL
#   for (xw in 1:4){
#   temp<-cbind(temp,matrix(l.neurons[[n]][[2]][xw] %*% x[,xw], ncol =1 ))
#   }
# l.neurons[[n]][[4]]<-temp[,1]+temp[,2]+temp[,3]+temp[,4]+l.neurons[[n]][[3]]
# }

for (generation in 
     
     