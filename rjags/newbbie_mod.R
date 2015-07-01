alpha=3
tau=2
y <- rnorm(10)
n <- length(y)


require(rjags)
#random effects over betaT (that is over the indicator variables)

data=list(y=y, n=n)
init=list(mean=rep(0,10))

modelstring="
  model {
    for (i in 1:n) {
      mean[i] ~ dbeta(alpha, b)
      y[i]~dnorm(mean[i],tau)
    }

    b ~ dgamma(c,d)
    d <- 1 
    c <- 0.01    
    alpha <- 3
    tau <- 0.01 
    
  }
"

model=jags.model(textConnection(modelstring),
                 data=data,inits=init)

update(model,n.iter=1000)
output=coda.samples(model=model,
                    variable.names=c("mean"),
                    n.iter=10000,thin=1)
print(summary(output))
plot(output)


