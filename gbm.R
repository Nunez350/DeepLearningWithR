setwd("~/Desktop/gbm-trainR-videos/data/")
?gbm

snp.data<-read.csv("~/msk-snp-cdgmp/cdg-data/cdgSNPmatrix-Jinyuan_roy.csv", sep =",", header = T, row.names = 1)
snp.by.gene <- split(snp.data, strtrim(rownames(snp.data), 8))
snp.by.gene.uniq <- lapply(snp.by.gene, function(x) unique(x))
snps<-lapply(snp.by.gene.uniq, function(x) t(x))
head(snps)
snps
cdgmp.data<-read.csv("cdgTable.csv2", sep =",", header =T)
m.cdgmp<-tapply(cdgmp.data$logcdg,cdgmp.data$strains,mean )
names(m.cdgmp)[16]<-gsub("pa_.+_", "", names(m.cdgmp)[16])
m.cdgmp<-m.cdgmp[rownames(snps[[1]])]
m.cdgmp
gene1<-snp.by.gene.uniq[[1]]
t <- as.matrix(ifelse(m.cdgmp < -1, 0,1 ))
tg<-t(gene1)
#library(Matrix)
mt<-tg#cbind(tg,t)
split<-sample(nrow(mt), floor(0.7*nrow(mt)))

train<-mt[split,]
test<-mt[-split,]

require(gbm)
require(dplyr)

test<-test[,-ncol(test)]
train<-train[,-ncol(train)]
end_trn = nrow(train)
all = rbind(train,test)
head(all)
dim(t)
dim(all)
t
gbm(~X1+X2+X3+X4+X5+X6,         # formula
    data=data
    }
end = nrow(all)
head(all)
colnames(all)
gsub(" ",",",noquote(split(colnames(all),","))," ",",")
?gsub
dim(all)
dim(all[1:end_trn,] )

ntrees = 5000
Model = gbm.fit( 
  x = all#[1:end_trn,] #dataframe of features
  , y = t #dependent variable
  #two ways to fit the model
  #use gbm.fit if you are going to specify x = and y = 
  #instead of using a formula
  #if there are lots of features, I think it's easier to specify 
  #x and y instead of using a formula
  
  
  , distribution = "bernoulli"
  #use bernoulli for binary outcomes
  #other values are "gaussian" for GBM regression 
  #or "adaboost"
  
  
  , n.trees = ntrees
  #Choose this value to be large, then we will prune the
  #tree after running the model
  
  
  , shrinkage = 0.01 
  #smaller values of shrinkage typically give slightly better performance
  #the cost is that the model takes longer to run for smaller values
  
  
  , interaction.depth = 3
  #use cross validation to choose interaction depth!!
  
  
  , n.minobsinnode = 10
  #n.minobsinnode has an important effect on overfitting!
  #decreasing this parameter increases the in-sample fit, 
  #but can result in overfitting
  
  , nTrain = round(end_trn * 0.8)
  #use this so that you can select the number of trees at the end
  
  # , var.monotone = c() 
  #can help with overfitting, will smooth bumpy curves
  
  , verbose = TRUE #print the preliminary output
)  
dim()
dim(t)





#look at the last model built
#Relative influence among the variables can be used in variable selection
summary(Model)
#If you see one variable that's much more important than all of the rest,
#that could be evidence of overfitting.

#optimal number of trees based upon CV
gbm.perf(Model)

#look at the effects of each variable, does it make sense?
?plot.gbm
for(i in 1:length(Model$var.names)){
  plot(Model, i.var = i
       , ntrees = gbm.perf(Model, plot.it = FALSE) #optimal number of trees
       , type = "response" #to get fitted probabilities
  )
}

###########################################









################ Make predictions ##################
#test set predictions
TestPredictions = predict(object = Model,newdata =all[(end_trn+1):end,]
                          , n.trees = gbm.perf(Model, plot.it = FALSE)
                          , type = "response") #to output a probability
#training set predictions
TrainPredictions = predict(object = Model,newdata =all[1:end_trn,]
                           , n.trees = gbm.perf(Model, plot.it = FALSE)
                           , type = "response")


#round the predictions to zero or one
#in general, don't do this!
#it was only because the answers in the comp had to be 0 or 1
TestPredictions = round(TestPredictions)
TrainPredictions = round(TrainPredictions)
#could also mess around with different cutoff values
#would need CV to determine the best


head(TrainPredictions, n = 20)
head(survived, n = 20)



#in sample classification accuracy
1 - sum(abs(survived - TrainPredictions)) / length(TrainPredictions) 
#depending upon the tuning parameters, 
#I've gotten this as high as 99%, but that model 
#resulted in lower test set scores


#to get predicted out of sample accuracy
#need to set aside a testing data set





#write the submission
submission = data.frame(PassengerId = 1:nrow(test), survived = TestPredictions)
write.csv(submission, file = "gbm submission.csv", row.names = FALSE)
#####################################################

N <- 1000

X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N) 
mu <- c(-1,0,1,2)[as.numeric(X3)]

SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)

# introduce some missing values
X1[sample(1:N,size=500)] <- NA
X4[sample(1:N,size=300)] <- NA
train<-tg
data <- data.frame(Y=t,X1=train[,1],X2=train[,2],X3=train[,3],X4=train[,4],X5=train[,5],X6= train[,6])
length(Y)
length(X1)
# fit initial model
gbm1 <-
  gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
      data=data,                   # dataset
      var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
      # +1: monotone increase,
      #  0: no monotone restrictions
      distribution="binomial",     # see the help for other choices
      n.trees=1000,                # number of trees
      shrinkage=0.05,              # shrinkage or learning rate,
      # 0.001 to 0.1 usually work
      interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
      bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
      train.fraction = 0.5,        # fraction of data for training,
      # first train.fraction*N used for training
      n.minobsinnode = 10,         # minimum total weight needed in each node
      cv.folds = 3,                # do 3-fold cross-validation
      keep.data=TRUE,              # keep a copy of the dataset with the object
      verbose=FALSE,               # don't print out progress
      n.cores=1)                   # use only a single core (detecting #cores is
# error-prone, so avoided here)

# check performance using an out-of-bag estimator
# OOB underestimates the optimal number of iterations
best.iter <- gbm.perf(gbm1,method="OOB")
print(best.iter)

