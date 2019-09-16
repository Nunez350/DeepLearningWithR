rm(list=ls(all=TRUE))
#single-layer multi neuron classifier
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

#list of neurons with activation function
l.neurons<-lapply(1:3, function(x) gneuron(x));names(l.neurons)<-c("z1","z2","z3")

historical.weight<-list();epoch<-NULL
neuron.upd.weights1<-data.frame();neuron.upd.weights2<-data.frame();neuron.upd.weights3<-data.frame()
for (epoch in 1:1000){
  #apply exponent function to all z
  z1e<-exp(unlist(l.neurons[[1]][1]));z2e<-exp(unlist(l.neurons[[2]][[1]]));z3e<-exp(unlist(l.neurons[[3]][[1]]))
  exp.act.list<-list(z1e,z2e,z3e)
  
  y1<-unlist(z1e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y2<-unlist(z2e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y3<-unlist(z3e)/(unlist(z1e)+unlist(z2e)+unlist(z3e))
  y=cbind(unlist(y1),unlist(y2),unlist(y3))
  
  e1=t1-y[,1]#calculates the error. Difference between the desired and predicted outpute1=t1-y[,1]
  e2=t2-y[,2]
  e3=t3-y[,3]
  e<-cbind(e1,e2,e3)
  
  updated.wx<-NULL
  gc<-NULL
  delta.w<-NULL
  n.r<-NULL
  sum.w<-NULL
  bias<-NULL
  delta.bias<-NULL
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
  
  r.bias<-NULL
  r.bias<-lapply(l.neurons, function(lx){
    matrix(unlist(lx[][3]))#nrow=1, ncol=4,byrow=T)
  })

  r.weights<-NULL
  r.weights<-lapply(l.neurons, function(lx){
  matrix(unlist(lx[][2]),nrow=1, ncol=4,byrow=T)
  
  })
  

    nwn<-list()
   for (n in 1:3){
     r.bias[[n]]<-cbind(r.bias[[n]],r.bias[[n]]-eta*delta.bias[[n]])
     
      temp<-data.frame(r.weights[[n]])
      for (ew in 1:4) { 
          for (i in 1:1){
          temp[i+1,ew]<-temp[i,ew]-eta*n.r[[n]][[ew]]
          }
        }
      nwn[[n]]<-temp
    }
    nwn
    r.bias[[1]]
    r.bias[[2]]
    r.bias[[3]]
    
    neuron.upd.weights1<-rbind(neuron.upd.weights1, nwn[[1]])
    neuron.upd.weights2<-rbind(neuron.upd.weights2, nwn[[2]])
    neuron.upd.weights3<-rbind(neuron.upd.weights3, nwn[[3]])
}

length(neuron.upd.weights1[,1])
length(neuron.upd.weights1[,2])
neuron.upd.weights1[,1]
plot(1:2000,neuron.upd.weights1[,1], ylim=c(-200,200), col="red", ylab= "weights", xlab="Epochs" )
lines(1:4000,neuron.upd.weights1[,2], col="blue")
lines(1:4000,neuron.upd.weights1[,3], col="green")
lines(1:4000,neuron.upd.weights1[,4], col="yellow")


lines(1:10000,neuron.upd.weights2, col="blue")
lines(1:10000,neuron.upd.weights3, col="blue")

#using weight and bias to predict species 
#r.bias[[1]], r.bias[[2]], r.bias[[3]]
r.bias
l.bias<- as.numeric(bias.list[10000])
fw1<-as.numeric(weights1[10000]); fw2<-as.numeric(weights2[10000])
y=1/(1+exp(-(fw1*x1 + fw2 * x2 + l.bias)))
round(y,1)
weights1[10000]
#target accuracy at each iteration
plot(1:10000,tacc)
ls()

str(nwn)[[]]
nwn<-data.frame(l.neurons[[1]][[2]],l.neurons[[2]][[2]],l.neurons[[3]][[2]])
gen<-data.frame()
for (g in 1:10){
  for (n in 1:3){
    r.bias[[n]]<-cbind(r.bias[[n]],r.bias[[n]]-eta*delta.bias[[n]])
    
    #temp<-data.frame(r.weights[[n]])
    for (ew in 1:4) { 
      rbind(nwn[[n]],do.call(data.frame, l.neurons[[n]][2])[[1]][[ew]]-eta*n.r[[n]][[ew]])
      #nwn[[n]]<-rbind(do.call(nwn[[n]],unlist(l.neurons[[n]][[2]])))
    }
    rbind(gen[[g]],do.call(data.frame,nwn))
    #rbind(nwn[g],nwn)
    #shaban[[n]]<-rbind(shaban,nwn)
  }
  #rbind(gen[[g]],do.call(data.frame,nwn))
  #rbind(nwn,do.call(data.frame,nwn))
}
dim(nwn)
do.call(data.frame,nwn)
do.call(data.frame, l.neurons[[1]][2])[[1]][[1]]
l.neurons[[n]][2][[1]]
l.neurons[[n]][2]

rbind(nwn[[n]],do.call(data.frame, l.neurons[[n]][2])[[1]])
 

nwn<-list(data.frame(l.neurons[[1]][[2]]),data.frame(l.neurons[[1]][[2]]),data.frame(l.neurons[[1]][[2]]))
nwn<-data.frame(l.neurons[[1]][[2]],l.neurons[[2]][[2]],l.neurons[[3]][[2]])
length(nwn[[2]])
#l.neurons[[1]][[2]]
temp<-data.frame()
for (n in 1:3){
  r.bias[[n]]<-cbind(r.bias[[n]],r.bias[[n]]-eta*delta.bias[[n]])

  #temp<-data.frame(r.weights[[n]])
  for (ew in 1:4) { 
   # for (i in 1:1){
   l.neurons[[n]][[2]][[ew]]=l.neurons[[n]][[2]][[ew]]-eta*n.r[[n]][[ew]]
  temp[[n]]=l.neurons[[n]][[2]][[ew]]
      #temp[i,ew]<-temp[i,ew]-eta*n.r[[n]][[ew]]
  #  }
  }
  
  nwn[[n]]<-rbind(do.call(nwn[[n]],unlist(l.neurons[[n]][[2]])))
}
nwn[[2]]
str(l.neurons[[n]][[2]])
rbind(nwn[[2]],do.call(data.frame, l.neurons[[1]][2])[[1]])

ff[,1]

n=1
do.call(data.frame, l.neurons[[1]][2])[[1]]
#nwn[[n]]<-
length(as.vector(nwn  ))
  rbind(nwn,as.vector(l.neurons[[n]][[2]]))

str(nwn[[n]])
str(l.neurons[[n]][[2]])
length
as.vector(l.neurons[[n]][[2]])
temp
class(l.neurons[[n]][[2]])
class(nwn)
nwn[[1]]
