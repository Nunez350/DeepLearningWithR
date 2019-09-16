#multi-layer neuron classifier
x<-iris
colnames(x)
t1<-c(rep(1,50),rep(0,50),rep(0,50))#first target
t2<-c(rep(0,50),rep(1,50),rep(0,50))#second target
t3<-c(rep(1,50),rep(0,50),rep(1,50))#third target
eta=0.1#the learning rate


weights1<-list();weights2<-list(); weights3<-list()
#x1<-iris[,1];x2<-iris[,2];x3<-iris[,3];x4<-iris[,4]
features<-list(x1<-iris[,1],x2<-iris[,2], x3<-iris[,3],x4<-iris[,4])

#generates random weights, initializes weights, random bias and activation
gneuron<-function(n){
  w<-runif(4, 1e-3, 1e-2)#random weights
  w1=w[1];w2=w[2];w3=w[3];w4=w[4]#initiating weights
  b<-runif(1)#random bias
  x1<-iris[,1];x2<-iris[,2]; x3<-iris[,3]; x4<-iris[,4]
  z=x1*w1 + x2*w2 + x3*w3 + x4*w4 +b
  h<-list(z,w[n],b)
  return(h)
}
gneuron(1)
#list of neurons with activation function
l.neurons<-lapply(1:3, function(x) gneuron(x))
names(l.neurons)<-c("z1","z2","z3")
l.neurons
z1e<-exp(unlist(l.neurons$z1[1]));z2e<-exp(unlist(l.neurons$z2[1]));z3e<-exp(unlist(l.neurons$z3[1]))
exp.act.list<-list(z1e,z2e,z3e)

y1<-unlist(z1e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y2<-unlist(z2e)/(unlist(z1e)+unlist(z2e)+unlist(z3e));y3<-unlist(z3e)/(unlist(z1e)+unlist(z2e)+unlist(z3e))
y=cbind(unlist(y1),unlist(y2),unlist(y3))
y[,1]
y
e1=t1-y[,1]#calculates the error. Difference between the desired and predicted outpute1=t1-y[,1]

e
#exp.list[1]/(exp.act.list[1]+exp.act.list[2]+exp.act.list[3])
#e=t1-unlist(y.list[1])

g1=lapply(unlist(features[[1]]),prod,-e)
g2=lapply(unlist(features[[2]]),prod,-e)
g3=lapply(unlist(features[[3]]),prod,-e)
g4=lapply(unlist(features[[4]]),prod,-e)
g.bias = -e

w1=w1 -eta * sum(g1)
# w2=w2 -eta * sum(g2)
# weights1[generation] <-w1
# weights2[generation] <-w2
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
