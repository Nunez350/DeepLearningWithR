rm(list=ls(all=TRUE))
#single-layer multi neuron classifier
x<-iris
x
colnames(x)
t1<-c(rep(1,50),rep(0,50),rep(0,50))#first target
t2<-c(rep(0,50),rep(1,50),rep(0,50))#second target
t3<-c(rep(0,50),rep(0,50),rep(1,50))#third target
eta=0.1#the learning rate

weights1<-list();weights2<-list(); weights3<-list(); weights4<-list()
#x1<-iris[,1];x2<-iris[,2];x3<-iris[,3];x4<-iris[,4]
#features<-list(x1<-iris[,1],x2<-iris[,2], x3<-iris[,3],x4<-iris[,4])

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
gneuron
#list of neurons with activation function
l.neurons<-lapply(1:3, function(x) gneuron(x))
names(l.neurons)<-c("z1","z2","z3")

#apply exponent function to all z
z1e<-exp(unlist(l.neurons[[1]][1]));z2e<-exp(unlist(l.neurons[[2]][[1]]));z3e<-exp(unlist(l.neurons[[3]][[1]]))
exp.act.list<-list(z1e,z2e,z3e)
exp.act.list
#computing y for each neuron
y1<-unlist(z1e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y2<-unlist(z2e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y3<-unlist(z3e)/(unlist(z1e)+unlist(z2e)+unlist(z3e))
y=cbind(unlist(y1),unlist(y2),unlist(y3))

e1=t1-y[,1]#calculates the error. Difference between the desired and predicted outpute1=t1-y[,1]
e2=t2-y[,2]
e3=t3-y[,3]

e<-cbind(e1,e2,e3)

#e<-as.data.frame(cbind(unlist(e1),unlist(e2),unlist(e3)))
#exp.list[1]/(exp.act.list[1]+exp.act.list[2]+exp.act.list[3])
#e=t1-unlist(y.list[1])

l.neurons[[1]]
l.neurons[[1]][[2]][[3]]
l.neurons[[1]][[3]]



X<-data.frame(unlist(l.neurons[[1]][[1]]), unlist(l.neurons[[2]][[1]]), unlist(l.neurons[[3]][[1]]))

df<-data.frame(x1=integer(), x2=integer(), x3=integer())      
updated.weights<-NULL
g<-NULL
sum.weights<-NULL
for (n in 1:3){
  for (r in 1:150){
    g<-rbind(g,prod(X[r,1], e[r,1]))
  }
  updated.weights<-cbind(updated.weights,g)
  g<-NULL
  sum.weights[n]<-sum(updated.weights[,n])
}

updated.weights
sum.weights

bias<-NULL
for(b in 1:3){
  
  bias<-append(bias,(-sum(e[,b])/150))
}
bias
plot(1:150, updated.weights[,1])





plot(1:10000,Dw ,ylim=c(-200,200), col="red", ylab= "weights", xlab="Epochs" )

# g1=lapply(unlist(features[[1]]),prod,-e)
# g2=lapply(unlist(features[[2]]),prod,-e)
# g3=lapply(unlist(features[[3]]),prod,-e)
# g4=lapply(unlist(features[[4]]),prod,-e)
# 
# w1=l.neurons[[2]][[2]][[1]] -eta * sum(unlist(g1))
# w2=l.neurons[[2]][[2]][[2]] -eta * sum(unlist(g2))
# w3=l.neurons[[2]][[2]][[3]] -eta * sum(unlist(g3))
# w4=l.neurons[[2]][[2]][[4]] -eta * sum(unlist(g4))
# 
# weights1[generation] <-w1
# weights2[generation] <-w2
# weights3[generation] <-w3
# weights4[generation] <-w4
# g.bias1 = -e1
# g.bias2 = -e2
# g.bias3 = -e3
# b = b - eta * sum(g.bias)
# bias.list[generation]<-b

# ncp<-length(which(abs(e) < 0.05))#number of predicted targets below erroe threshold
# tac<-ncp/100 #predicted target accuracy rate
# tacc[generation]<-tac
# }
# output<- data.frame(weights1=unlist(weights1), weights2=unlist(weights2),errors=unlist(bias.list))
# 
# #plotting changes of weights over epoch
# plot(1:1000,weights1, ylim=c(-200,200), col="red", ylab= "weights", xlab="Epochs" )
# lines(1:1000,weights2, col="blue")
# 
# #using weight and bias to predict species 
# l.bias<- as.numeric(bias.list[1000])
# fw1<-as.numeric(weights1[1000]); fw2<-as.numeric(weights2[1000])
# y=1/(1+exp(-(fw1*x1 + fw2 * x2 + l.bias)))
# round(y,1)
# 
# #target accuracy at each iteration
# plot(1:1000,tacc)
# 
# #classifing two species
# class.species<-function(){
#   sp0<-x[which(x$Species == "versicolor"),]
#   sp1<-x[which(x$Species == "virginica"),]
#   plot(sp0$Sepal.Length, sp0$Petal.Length, col="red", ylim = c(2.9,6.4), ylab = "Petal length", xlab="Sepal length")
#   points(sp1$Sepal.Length, sp1$Petal.Length, col="blue")
#   legend("topleft", legend = c("versicolor","virginica"), col=c("red","blue"),pch=1, horiz = F, cex = .8) 
#   slope=-(w1/w2);intercept=-b/w2
#   abline(slope,intercept)
# }
# class.species()
cat<-0;Dw<-list(); Dbias<-list()

dim(l.neurons)
sum(x[,1])
x[,1]
e
Dw
l.neurons[[1]][[1]]     
df<-data.frame(x1=integer(), x2=integer(), x3=integer())      
df<-lapply(1:3, function(dx){ 
  temp<-unlist(l.neurons[[1]])
  
  df<-cbind(df,temp)
  temp<-NULL
  
  
}
)
df<-lapply(1:3, function(dx) unlist(l.neurons[[1]]))
df

l.neurons
dx<-data.frame()
dx
for (i in 1:3){
  dx<-cbind(dx, unlist(l.neurons[[i]][[1]]))
}
#Z<-data.frame(unlist(l.neurons[[1]][[1]]), unlist(l.neurons[[2]][[1]]), unlist(l.neurons[[3]][[1]]))
#Z

#df<-data.frame(x1=integer(), x2=integer(), x3=integer())  