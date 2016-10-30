### Load packages
library(lars)
library(glmnet)
library(MASS)

genNorm <- function(x, p) {
  sum(abs(x)^p)^(1/p)
  
}

############################################
# Create synthetic data
############################################
# Set parameters 
nbFeatures<-100
nbSamples<-200
nbRelevantIndices<-3
snr <- 1      # signal noise ration
nbLambdas<-10
nlambda<-nbLambdas

sigma<-matrix(0,nbFeatures,nbFeatures)
for (i in (1:nbFeatures))
  for (j in (1:nbFeatures))
    sigma[i,j]<-0.99^(abs(i-j))

X <- mvrnorm(n = nbSamples, mu=rep(0, nbFeatures), Sigma=sigma, tol = 1e-6, empirical = FALSE)
relevantIndices <- sample(1:nbFeatures, nbRelevantIndices, replace=F)

X <- scale(X) #bug fixed, replaced 'scale(X)' by 'X<-scale(X)'
beta <- rep(0, times=nbFeatures)
beta[relevantIndices] <- runif(nbRelevantIndices,0,1)

epsilon<-matrix(rnorm(nbSamples), nbSamples, 1)

simulatedSnr<-(norm(X %*% beta,type="F"))/(norm(epsilon,type="F"))   # simulated Signal Noise Ratio
epsilon <- epsilon*(simulatedSnr/snr)  
 
Y <- X %*% beta + epsilon;
Y <- scale(Y) #bug fixed, replaced 'scale(Y)' by 'Y<-scale(Y)'
#hist(Y)   # check type of response to set the family in glmnet


# ### launch lars 
# la <- lars(X,Y,intercept=TRUE, max.steps=1000, use.Gram=FALSE, trace=T, normalize=F)
# coLars <- as.matrix(coef(la,,mode="lambda"))
# #print(round(coLars,4))
# #plot(la)
# coeffs = 0
# lastIter <- dim(coLars)[1]
# coeffs <- rbind(coLars[lastIter, ], beta)
# barplot(coeffs, col=c("blue", "red"), leg=c("predicted", "true"))


### Launch glmnet with lambda=1/2*lambda_lars 
# penalty = rep(1, nbFeatures)
# lambda=0.5*la$lambda
# glm2 <- glmnet(X,Y,family="gaussian", thresh=1e-16, 
#                standardize=F, type.gaussian="covariance", 
#                penalty.factor=penalty)

glm2 <- cv.glmnet(X,Y, standardize=F, type.gaussian="covariance")
coGlm2 <- as.matrix(t(coef(glm2, mode="lambda")))
#print(round(coGlm2,4))
#plot(glm2)

testSet <- sample(1:nbSamples, replace=F, size=10)
predY <- predict(glm2, newx=X[testSet, ], s="lambda.min")


coeffs = 0 
lastIter <- dim(coGlm2)[1]
coeffs <- rbind(coGlm2[lastIter, -1], beta)
barplot(coeffs, col=c("blue", "red"), leg=c("predicted", "true"))

cor(coeffs[1,], coeffs[2,])

# coefficient error
coeffError <- genNorm(coeffs[2,]-coeffs[1,],1)

# estimation error
estError <- genNorm(predY-Y[testSet], 2) 

cor(predY, Y[testSet])





####################################################################################
## visualise differences in LARS and GLMNET coefficients 
coeffs <- rbind(round(coLars[dim(coLars)[1], ],4) , round(coGlm2[dim(coGlm2)[1],-1],4 ))
barplot(coeffs, col=c("red", "green"), 
        leg=c("lars","glmnet"), 
        main="Glmnet vs Lars", 
        ylab="|beta|", xlab="Coefficients")

#cor(coGlm2[13, -1], coLars[14, ])


##########################
#edit: The rest of the program is not relevant any more, as the issue of the intercept was due to a bug.    
##########################

### Remove intercept computed by glmnet for some lambda (adjusted by hand, depending on output) 
###Y<-Y--0.0863
###Y<-Y/norm(Y,"F")

### Launch Lars 
###la <- lars(X,Y,intercept=TRUE, max.steps=1000, use.Gram=FALSE)
###coLars <- as.matrix(coef(la,,mode="lambda"))
###print(round(coLars,4))

### Launch Glmnet with lambda=1/2*lambda_lars 
###glm2 <- glmnet(X,Y,family="gaussian",lambda=0.5*la$lambda,thresh=1e-16)
###coGlm2 <- as.matrix(t(coef(glm2,mode="lambda")))
###print(round(coGlm2,4))