# simulate SNP-cdg correlation
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



r <- 0.8 # desired correlation coefficient
sigma <- matrix(c(1,r,r,1), ncol=2) # var-covariance matrix
?chol
s <- chol(sigma) # choleski decomposition
n <- 100 # number of random deviates (data points)
z <- s %*% matrix(rnorm(n*2), nrow=2) # 100 correlated normally distributed deviates with cor(x,y)=r
u <- pnorm(z) # get probabilities for each deviates
snp.states <- qbinom(u[1,], 1, 0.5) # discretize the 1st vector of probabilities into 0/1 with Bernoulli trial
idx.0 <- which(snp.states == 0); # indices for "0"
idx.1 <- which(snp.states == 1); # indices for "1"
# boxplots with stripcharts
boxplot(u[2,] ~ snp.states, main="cor=0.8", xlab="SNP states", ylab="CDG level")
stripchart(u[2,] ~ snp.states, vertical=T, pch=1, method="jitter", col=2,  add=T)


r <- 0.5 # desired correlation coefficient
sigma <- matrix(c(1,r,r,1), ncol=2) # var-covariance matrix
s <- chol(sigma) # choleski decomposition
s
n <- 100 # number of random deviates (data points)
z <- s %*% matrix(rnorm(n*2), nrow=2) # 100 correlated normally distributed deviates with cor(x,y)=r
z
u <- pnorm(z) # get probabilities for each deviates
u
snp.states <- qbinom(u[1,], 1, 0.5) # discretize the 1st vector of probabilities into 0/1 with Bernoulli trial
?qbinom
snp.states
idx.0 <- which(snp.states == 0); # indices for "0"
idx.1 <- which(snp.states == 1); # indices for "1"
# boxplots with stripcharts
u[2,]
boxplot(u[2,] ~ snp.states, main="cor=0.5", xlab="SNP states", ylab="CDG level")
stripchart(u[2,] ~ snp.states, vertical=T, pch=1, method="jitter", col=2,  add=T)




library(SimPhe)

x1 <- rnorm(4000, mean = 5, sd = 10)
x2 <- rnorm(4000, mean = 10, sd = 30)
x <- matrix(cbind(x1, x2), ncol = 2)
# test original correlation
cor.test(x[, 1], x[, 2])
# correlation matrix
corM <- matrix(c(1, 0.6, 0.6, 1), ncol = 2)
# standard deviation matrix
sdM <- matrix(c(10, 0, 0, 30), ncol = 2)
# build correlation
x.new <- build.cor.phe(x, corM, sdM)
# check mean and standard deviation of new data set
apply(x.new, 2, mean)
apply(x.new, 2, sd)
# test correlation
cor.test(x.new[, 1], x.new[, 2])
x.new[, 1]

x1 <- rnorm(4000, mean = 5, sd = 10)
x2 <- rnorm(4000, mean = 10, sd = 30)
x <- matrix(cbind(x1, x2), ncol = 2)
build.sd.matrix(x)





library(scrime)
simulateSNPglm()
sim.snp<-simulateSNPglm(30, 50, list.ia = NULL, list.snp = NULL, 
               beta0 = -0.5, beta = 1.5, maf = 0.25, sample.y = TRUE, p.cutoff = 0.5, rand = NA)
sim.snp
