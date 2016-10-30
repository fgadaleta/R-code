require(ElemStatLearn)
#data(galaxy)
#y = galaxy[,5]
#x = as.matrix(galaxy[,-5])
set.seed(13)

simbeta = function(nbeta, main) {
  truebeta = rnorm(nbeta,0,sd=0.01)
  eff = sample(1:10, main, replace=F)
  truebeta[eff] = rnorm(main, 0.7, 0.01)
  return(truebeta)
}

# synthetic data w/wo correlations ######
n=4;p=20
x=matrix(rnorm(n*p),ncol=p)

# generate correlated covariates x
grp = 3; vpg = 4
corvars = sample(1:p, size=grp*vpg, replace=F)
for (i in 1:grp){
  from = (i-1)*vpg+1
  to = (i*vpg)
  g1 = corvars[from:to]
  x[, g1] = x[, g1[1]] + rnorm(1, mean=0,sd=1)
}

heatmap(cor(x))

#truebeta = simbeta(p, 3)
truebeta = rnorm(p, 0, 0.001)
truebeta[11] = rnorm(1, 1, 0.01)
truebeta[13] = rnorm(1, 1, 0.01)
truebeta[16] = rnorm(1, 1, 0.01)
tbeta = round(truebeta, digits = 2)

alpha=3;tau=2
eps=rnorm(n,0,1/sqrt(tau))
y=alpha+as.vector(x%*%tbeta + eps)

mod=lm(y~x)
print(summary(mod))
#########################################


################################

which(tbeta>0.9)

# linear regression
linmod = lm(y~x)
plot(linmod)

# ordinary least square
x0 = cbind(rep(1, length(y)), x)
ols = lm(y~0+x0)
plot(ols)



# t(x)*x*B = t(x)y
solve(t(x)%*%x, t(x)%*%y)

#scale data 
y_scaled = y-mean(y)
W = sweep(x,2,colMeans(x))
solve(t(W)%*%W, t(W)%*%y)    # OLS

# ridge regression with lambda = 100 
lambda = 100
rrpred = solve(t(W) %*% W + lambda*diag(dim(W)[2]), t(W)%*%y_scaled)
sort(round(rrpred, 2), decreasing=T, index=T)$ix


###################################################################
###### define a loss function and optimise directly        ########
###################################################################
loss_en<-function(beta) {
   eps=y_scaled-W%*%beta
   return(sum(eps*eps)+lambda1*sum(abs(beta))+lambda2*sum(beta*beta)) 
}

loss_lasso<-function(beta) {
  eps = y_scaled - W%*%beta
  return(sum(eps*eps) + lambda1*sum(abs(beta))) 
}

loss_ridge<-function(beta) {
  eps = y_scaled - W%*%beta
  return(sum(eps*eps) + lambda2*sum(beta*beta)) 
}

loss_fused<-function(beta) {
  eps = y_scaled - W%*%beta  # loss 
  
  p = length(beta)
  pen = {}
  for(i in 2:p){ 
    tmp = abs(beta[i]-beta[i-1])
    pen = c(tmp, pen)
  }

  return(sum(eps*eps) + lambda2*sum(abs(beta)) + lambda1*sum(pen)) 
}




lambda1=2    # lasso penalty  
lambda2=1    # ridge regression penalty  

betaen = optim(rep(0,dim(W)[2]),loss_en,control=list(maxit=10000,reltol=1e-12))$par
which(round(betaen, 2) != 0)

betalasso = optim(rep(0,dim(W)[2]),loss_lasso,control=list(maxit=10000,reltol=1e-12))$par
which(round(betalasso, 2) != 0)

betaridge = optim(rep(0,dim(W)[2]),loss_ridge,control=list(maxit=10000,reltol=1e-12))$par
which(round(betaridge, 2) != 0)

betafused = optim(numeric(dim(W)[2]),loss_fused,control=list(maxit=10000,reltol=1e-12))$par
which(round(betafused, 2) != 0)
plot(betafused)

tmp = fusedlasso1d(y = y_scaled, X = W, maxsteps = 10000)
plot(tmp$beta[,3])
tmp$lambda
tmp$
#round(betapred, digits=2)
#sort(round(abs(betapred), digits=2), decreasing=T, index=T)$ix

###################################################################
###  use enet package only for loss functions of this form      ###
###################################################################
require(elasticnet)
obj = enet(W,y_scaled, lambda=lambda2, normalize=F)
coefs = predict(obj, s=lambda1, mode="penalty",type="coefficients")$coefficients
which(coefs!=0)
#sort(round(coefs, digits=2), decreasing=T, index=T)$ix[1:10] 

###################################################################
###  use glmnet                                                 ###
###################################################################
library(glmnet)
penalties = rep(1,times = p)        # if 0, no shrinkage (var always included) (default is 1)
penalties[9]= penalties[13] = 1000  # if high, exclude these ones 

fit = glmnet(W, y_scaled, standardize=F, intercept=F, lambda = 0.8)
fit = glmnet(W, y_scaled, standardize=F, intercept=F, lambda = 0.8, penalty.factor = penalties) 
which(fit$beta != 0)

#cvfit = cv.glmnet(x = W, y = y_scaled)
#fit = glmnet(W, y_scaled, standardize=F, intercept=F, lambda = cvfit$lambda.min)
#coefs = coef(fit)
#which(coefs!=0)

plot(fit)
cormat = cor(W)
sort(cormat[26,], decreasing=T, index=T)





library(grpreg)
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="grLasso")
fit <- grpreg(x,y, intercept=F)   # group set up automatically
fit$group
plot(fit)
betas = fit$beta
select(fit, "AIC")

AIC(fit)
# TODO group variables by correlation

library(grplasso)

data <- read.delim("../data/gnw_50_3/gnw_50_3-1_nonoise_multifactorial.tsv")
p<-ncol(data)
n<-nrow(data)
data.scale <- scale(data)
names(data.scale) = names(data)
z=1
noti <- (1:p)[-z]
yi <- data.scale[, z]
xi <- data.scale[, noti]    
pbeta = {}     # vector of predicted beta 
fit = cv.glmnet(x=xi, y=yi, standardize=F, intercept=F, alpha=1, type.measure="mae")  # mse    
cvlambda = fit$lambda.min 
fit = glmnet(xi, yi, standardize=F, intercept=F, alpha=1, lambda=cvlambda) 
singlecoeffs <- as.matrix(t(coef(fit, s="lambda.min")))
singlecoeffs = singlecoeffs[,-1]                      # remove intercept
maincoeff = sort(abs(singlecoeffs), decreasing=T, index.return=T)$ix
effects = singlecoeffs[maincoeff]
mainlasso  = which(effects != 0 )


corrmat <- cor(xi)
group <- rep(1,times=p-1)
### TODO set groups by correlation 
#(highly correlated variables belong to same group)
current = corrmat[1,1]  
for(i in 1:(p-1)) {
  g = which(corrmat[i,] > 0.7)
  g = g[-i]
  group[i] = i
}

#clustering on correlation might help 



group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8,1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8,
           1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8,8)

fit = cv.grpreg(xi, yi, group=group)
cvlambda = fit$lambda.min
fitgrp <- grpreg(xi,yi, lambda=cvlambda, group=group)
betas = as.matrix(fitgrp$beta)
betas = betas[-1,]
mainbetas = sort(abs(betas), decreasing=T, index.return=T)$ix
maingrlasso = which(betas != 0)


library(hierNet)
fitcv = hierNet.path(xi,yi)
singlefit = hierNet.cv(fitcv, xi,yi)
singlefit = hierNet(x= xi, y= yi, strong=T, lam=singlefit$lamhat, niter=100, maxiter=1000)

#TODO main effects 
#singlefit$bp # positive main effect
#singlefit$bn # negative main effect
singlecoeffs = singlefit$bp - singlefit$bn # overall main effect estimated coefficients
singleinter = round(singlefit$th,2)
heatmap(singleinter)
hierlasso = round(singlecoeffs,2)

library(ggplot2)
qplot(hierlasso, geom="density") + theme_bw()
which(hierlasso>0)






