#!/usr/bin/Rscript
library(lassoshooting)
# make data 
rm(list = ls(all = TRUE)) # make sure previous work is clear
ls()
#x0 <- c(1,1,1,1,1) # column of 1's
#x1 <- c(1,2,3,4,5) # original x-values
# create the x- matrix of explanatory variables
#x <- as.matrix(cbind(x0,x1))
p = 70; n = 300
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
# create the y-matrix of dependent variables
#y <- as.matrix(c(3,7,5,11,14))
y <- matrix(rnorm(n), nrow=n)
#y <- as.matrix(rnorm())
m <- nrow(y)
# implement feature scaling
#x.scaled <- x
#x.scaled[,2] <- (x[,2] - mean(x[,2]))/sd(x[,2])

# analytical results with matrix algebra
analytical <- solve(t(x)%*%x)%*%t(x)%*%y # w/o feature scaling
#solve(t(x.scaled)%*%x.scaled)%*%t(x.scaled)%*%y # w/ feature scaling
# results using canned lm function match results above
#summary(lm(y ~ x[, 2])) # w/o feature scaling
#summary(lm(y ~ x.scaled[, 2])) # w/feature scaling

# define the gradient function dJ/dtheata: 1/m * (h(x)-y))*x where h(x) = x*theta
# in matrix form this is as follows:
grad <- function(x, y, theta) {
  gradient <- (1/m)* (t(x) %*% ((x %*% t(theta)) - y))
  return(t(gradient))
}


# define gradient descent update algorithm
grad.descent <- function(x, maxit){
  thetas = matrix(0,ncol=dim(x)[2], nrow=maxit)
  theta <- matrix(0, nrow=1, ncol=dim(x)[2])
  numcoord = dim(x)[2]
  
  alpha = .1 # set learning rate
  
  for (i in 1:maxit) {
    theta <- theta - alpha*grad(x, y, theta)   
    thetas[i, ] = theta 
  }
  return(thetas)
}

coord.descent <- function(x, maxit, lambda=0.1){
  p = dim(x)[2]
  theta <- matrix(0, ncol=p, nrow=maxit) # Initialize the parameters
  
  for(k in 1:maxit){
    for(i in 1:p) {
      iter = ifelse(k>1, yes=k-1, 1)
      theta[k,i] = t(x[,i])%*% (y - x[,-i]%*%theta[iter,-i])/ (t(x[,i])%*%x[,i])
      theta[k,i] = softthresh(theta[k,i], lambda)
    }
  }
  return(theta)
}

maxiter = 20
out_gd = grad.descent(x, maxiter)
out_cd = coord.descent(x, maxiter, lambda=0)


library(ggplot2)
out1 = data.frame(iter=1:maxiter ,p2 = out_cd[,2])
out2 = data.frame(iter=1:maxiter ,p2 = out_gd[,2])
anal = data.frame(iter=1:maxiter, p2 = analytical[2,])

png(filename="./gd_vs_cd.png", width=1200, height=800)

ggplot(out1,aes(iter,p2)) + geom_line(aes(color="Coordinate descent"), size=1.5) +
  geom_line(data=out2,aes(color="Gradient descent", size=1)) + labs(color="") +
  geom_line(data=anal, aes(color="Analytical", size=1))

dev.off()

# ----------------------------------------------------------------------- 
# cost and convergence intuition
# -----------------------------------------------------------------------
# typically we would iterate the algorithm above until the 
# change in the cost function (as a result of the updated b0 and b1 values)
# was extremely small value 'c'. C would be referred to as the set 'convergence'
# criteria. If C is not met after a given # of iterations, you can increase the
# iterations or change the learning rate 'alpha' to speed up convergence

# get results from gradient descent
# beta <- grad.descent(x,1000)
# beta <- coord.descent(x, 1000)
# 
# # define the 'hypothesis function'
# h <- function(xi,b0,b1) {
#   b0 + b1 * xi 
# }
# 
# # define the cost function   
# cost <- t(mat.or.vec(1,m))
# for(i in 1:m) {
#   cost[i,1] <-  (1 /(2*m)) * (h(x[i,2], beta[1000,1], beta[1000,2]) - y[i,])^2 
# }
# 
# totalCost <- colSums(cost)
# print(totalCost)
# 
# # save this as Cost1000
# cost1000 <- totalCost
# 
# # change iterations to 1001 and compute cost1001
# beta <- (grad.descent(x,1500))
# cost <- t(mat.or.vec(1,m))
# for(i in 1:m) {
#   cost[i,1] <-  (1 /(2*m)) * (h(x[i,2],beta[1,1],beta[1,2])- y[i,])^2 
# }
# cost1001 <- colSums(cost)
# 
# # does this difference meet your convergence criteria? 
# print(cost1000 - cost1001)
# 


