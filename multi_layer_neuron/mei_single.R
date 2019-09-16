cdg <- read.csv("cdg-data/cdgTable.csv2", header = T, na.strings = T)
gene1 <- read.csv("cdg-data/gene-1.snp.csv2", header = T, row.names = 1)


snps <- t(gene1)
snps[,1]
cdg.mean <- tapply(cdg$logcdg, cdg$strains, mean)
t <- as.matrix(ifelse(cdg.mean[row.names(snps)] < -1, 0,1 ))
nrow(t)

snps <- cbind(t, snps)
row.names(snps)

w <- runif(32, 1e-3, 1e-2)
b <- runif(1)

eta <- 0.1

weights <- matrix(0,nrow = 1000, ncol = 32)
accuracy <- numeric()
snps[,1]
dim(snps)
for(epoch in 1:1000){
  a <- snps[,2:33] %*% w
  y <- 1/(1+exp(-a-b))
  e <- snps[,1] - y
  w <- w - eta * -colSums(snps[,2:33] * e[,1])
  b <- b - eta * sum(-e)
  for(k in 1:ncol(weights)){
    weights[epoch,k] <- w[k]
  }
  accuracy <- c(accuracy, length(which(round(y) == t))/30)
}

plot(accuracy, type = "p")

