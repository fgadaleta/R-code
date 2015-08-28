# Generate data 

m=50   # columns or features 
n=100  # rows or observations
frac.gaps <- 0.5 # the fraction of data with NaNs
N.S.ratio <- 0.25 # the Noise to Signal ratio for adding noise to data

x <- (seq(m)*2*pi)/m
t <- (seq(n)*2*pi)/n

#True field
Xt <- 
  outer(sin(x), sin(t)) + 
  outer(sin(2.1*x), sin(2.1*t)) + 
  outer(sin(3.1*x), sin(3.1*t)) +
  outer(tanh(x), cos(t)) + 
  outer(tanh(2*x), cos(2.1*t)) + 
  outer(tanh(4*x), cos(0.1*t)) + 
  outer(tanh(2.4*x), cos(1.1*t)) + 
  tanh(outer(x, t, FUN="+")) + 
  tanh(outer(x, 2*t, FUN="+"))

Xt <- t(Xt)   # data matrix nxm


# PCA 
res <- prcomp(Xt, center = TRUE, scale = FALSE)
names(res)

#plot(cumsum(res$sdev^2/sum(res$sdev^2))) #cumulative explained variance

##################################################################
# reconstruction of original data from lower number of features  #
##################################################################
pc.use <- 3                        # num of principal components to use 
Xt_projected <- res$x[,1:pc.use]   # projection of original data onto PCs
Xt_reconstructed <- Xt_projected %*% t(res$rotation[,1:pc.use])

# add the center (and re-scale) back to original data
if(res$scale != FALSE){
  Xt_reconstructed <- scale(Xt_reconstructed, center = FALSE , scale=1/res$scale)
}
if(res$center != FALSE){
  Xt_reconstructed <- scale(Xt_reconstructed, center = -1 * res$center, scale=FALSE)
}
dim(Xt_reconstructed); dim(Xt)

# how much space we save
# total size = size of projected data + size of rotation matrix 
compressionRatio = as.numeric(object.size(Xt)/(object.size(res$x[,1:pc.use]) + object.size(res$rotation[,1:pc.use])))


######################################
## visualization                    ##
######################################
RAN <- range(cbind(Xt, Xt_reconstructed))
BREAKS <- seq(RAN[1], RAN[2],,100)
COLS <- rainbow(length(BREAKS)-1)
par(mfcol=c(1,2), mar=c(1,1,2,1))
image(Xt,xlab="", ylab="", xaxt="n", yaxt="n", breaks=BREAKS, col=COLS)
box()
image(Xt_reconstructed, xlab="", ylab="", xaxt="n", yaxt="n", breaks=BREAKS, col=COLS)
box()
rm(RAN)
rm(BREAKS)
rm(COLS)






######################################
#### DIY version                  ####
######################################
#Z <- as.matrix(scale(Xt, center=TRUE, scale=FALSE))

# compute mean per column
mean <- colMeans(Xt)      
#Xt_center <- scale(Xt, center = 1 * mean, scale=FALSE) 

# make matrix of means to be subtracted to original data
ones = rep(1, nrow(Xt))   
x_mean = ones %*% t(colMeans(Xt))

# original data get centered
Xt_center <- Xt - x_mean

# eigen values/vectors of correlation matrix 
K <- eigen(cor(Xt_center))
numpc = 5

# projection of original data onto the principal components
P <- Xt_center %*% -K$vec[,1:numpc]  
# reconstruction from projected data back 
Z_hat <- P %*% t(K$vectors[,1:numpc])

# add mean back 
Z_hat <- Z_hat + x_mean
#Z_hat <- scale(Z_hat, center = 1 * mean, scale=FALSE)   # more space efficient no need to store a matrix of means

# visualize comparison original/reconstructed 
RAN <- range(cbind(Xt, Z_hat))
BREAKS <- seq(RAN[1], RAN[2],,100)
COLS <- rainbow(length(BREAKS)-1)
par(mfcol=c(1,2), mar=c(1,1,2,1))
image(Xt, xlab="", ylab="", xaxt="n", yaxt="n", breaks=BREAKS, col=COLS)
box()
image(Z_hat, xlab="", ylab="", xaxt="n", yaxt="n", breaks=BREAKS, col=COLS)
box()
rm(RAN)
rm(BREAKS)
rm(COLS)