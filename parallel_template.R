require(snowfall)


p = 2000; n = 200
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
#df <- data.frame(x)

# split data to parallelise
mats = {}
for(i in 1:p){
  yi = x[,i]
  xi = x[,-i]
  mi = cbind(yi, xi)
  mats[[i]] = mi
}


singlefit <- function(dat) {
  fit = glmnet(y=dat[,1], x=dat[,-1], alpha=1, intercept=F, lambda=0.1)
  singlecoeffs <- as.matrix(t(coef(fit)))
  singlecoeffs <- singlecoeffs[,-1]
  return(singlecoeffs)
}

start = Sys.time()

sfInit(parallel = T, cpus = 4, type="SOCK")
sfLibrary(glmnet)
res <- sfLapply(mats, singlefit)
sfStop()

end = Sys.time()
elapsed = end-start
elapsed*60

