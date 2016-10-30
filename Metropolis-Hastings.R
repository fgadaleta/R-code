f <- function(mu, alpha, beta, gamma, delta, tau, y) {
  c<-log( beta^(gamma-1)*exp(-delta*beta) )
  b <- 0
  
  for(i in 1:length(y)) {
    b <- b -(y[i]-mu[i])^2 + log( mu[i]^(alpha-1)*(1-mu[i])^(beta-1)) 
  }
  
  a <- b+c
  return(a)
}



mu <- rnorm(10,0)
y  <- rnorm(10, 1)
alpha <- 3
beta <- 3
gamma <- 2
delta <- 1
tau <- 0.01



f(mu, alpha, beta, gamma, delta, tau, y)


mu.draws <- numeric(5000)
mu.last <- 1

for(i in 1:5000) {
  mu.cand <- runif(1,0,1)
  ratio <- f(mu.cand, alpha, beta, gamma, delta, tau, y)/ f(mu.last, alpha, beta, gamma, delta, tau, y)
  if(runif(1)< min(ratio,1)) {
    mu.last <- mu.cand
  }
  
  mu.draws[i] <- mu.last
  
  }


