co# simulate SNP-cdg correlation
r <- 0.2 # desired correlation coefficient
sigma <- matrix(c(1,r,r,1), ncol=2) # var-covariance matrix
s <- chol(sigma) # choleski decomposition
n <- 100 # number of random deviates (data points)
z <- s %*% matrix(rnorm(n*2), nrow=2) # 100 correlated normally distributed deviates with cor(x,y)=r
u <- pnorm(z) # get probabilities for each deviates
snp.states <- qbinom(u[1,], 1, 0.5) # discretize the 1st vector of probabilities into 0/1 with Bernoulli trial
idx.0 <- which(snp.states == 0); # indices for "0"
idx.1 <- which(snp.states == 1); # indices for "1"
# boxplots with stripcharts
boxplot(u[2,] ~ snp.states, main="cor=0.2", xlab="SNP states", ylab="CDG level")
stripchart(u[2,] ~ snp.states, vertical=T, pch=1, method="jitter", col=2,  add=T)
u
simulated.training<-u[2,]
snp.states


t1 <- as.matrix(ifelse(simulated.training > -.5, 1,0 )) # high
t2 <- as.matrix(ifelse(simulated.training < -.5 & simulated.training > -1.5, 1,0)) # medium 
t3 <- as.matrix(ifelse(simulated.training  < -1.5 , 1,0 )) #low
targets<-cbind(t1,t2,t3)
targets


snps
