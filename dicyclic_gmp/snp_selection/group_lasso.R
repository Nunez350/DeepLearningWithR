library("grplasso")
data(splice)
str(splice)
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]
fit <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 10,
                contrasts = contr, standardize = TRUE)
splice$contr
fit <- grplasso(t~ ., data = as.data.frame(t(snp.by.gene.uniq)), model = LogReg(), lambda = 10,contrasts = contr, standardize = TRUE)
snp<-as.data.frame(t(snp.by.gene.uniq))
snp<-as.matrix(t(snp.by.gene.uniq))
t<-as.matrix(m.cdgmp)


fit <- grplasso(t~ ., data = man), model = LogReg(), lambda = 10,contrasts = contr, standardize = TRUE)
snp<-cbind(snp,t)
snp
as.data.frame(t(snp.by.gene.uniq))
cbind(t(snp.by.gene.uniq)[[1]],t(snp.by.gene.uniq)[[1]])
cbind(t(snp.by.gene.uniq[[2]]),t(snp.by.gene.uniq[[2]]))
man<-NULL
for (i in 1:50){
  
  man<-cbind(man,t(snp.by.gene.uniq[[i]]))
}
t
#man<-
  man[,1078]

dim(man)
#man<-t(snp.by.gene.uniq[[2]])
man<-cbind(man,t)
man[,t]
man<-as.data.frame(man)
grplasso(as.matrix(t(snp.by.gene.uniq[[1]])),as.matrix(t))
snp.by.gene.uniq
as.matrix(t)
as.matrix(t)
as.matrix(man)
str(man)
man[t]








## Use the Logistic Group Lasso on the splice data set
data(splice)

## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model 
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20,
                       contrasts = contr, center = TRUE, standardize = TRUE)

## Perform the Logistic Group Lasso on a random dataset
set.seed(79)

n <- 50  ## observations
p <- 4   ## variables

## First variable (intercept) not penalized, two groups having 2 degrees
## of freedom each

index <- c(NA, 2, 2, 3, 3)

## Create a random design matrix, including the intercept (first column)
x <- cbind(1, matrix(rnorm(p * n), nrow = n))
colnames(x) <- c("Intercept", paste("X", 1:4, sep = ""))

par <- c(0, 2.1, -1.8, 0, 0)
prob <- 1 / (1 + exp(-x %*% par))
mean(pmin(prob, 1 - prob)) ## Bayes risk
y <- rbinom(n, size = 1, prob = prob) ## binary response vector

## Use a multiplicative grid for the penalty parameter lambda, starting
## at the maximal lambda value
lambda <- lambdamax(x, y = y, index = index, penscale = sqrt,
                    model = LogReg()) * 0.5^(0:5)

## Fit the solution path on the lambda grid
fit2 <- grplasso(mat, y = groups, index = index, lambda = lambda, model = LogReg(),
                penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))
summary(fit2)
grplasso(x, y = y, index = index)
## Plot coefficient paths
plot(fit)
index
length(groups)
dim(x)
length(y)
length(index)
dim(mat)
length(groups)
index=c(1:7)
groups<-y[1:7]
mat<-as.matrix(cbind(rep(1,30),t(snp.by.gene.uniq[[1]])))
mat
dim(mat)
grplasso(mat,groups, index)
x
dim(mat)
length(groups)
length(index)       
length(groups)
length(index)
dim(mat)
ncol(mat)
mat
groups <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))
grplasso(as.matrix(cbind(rep(1,30))),t(snp.by.gene.uniq[[1]]),as.matrix(t))
str(y)
groups<-y[1:7]
ncol(x)
index=colnames(snp.by.gene.uniq[[1]])
index=c(1:30)
