# hidden markov model of a dishonest casino
# A transition matrix
# E emission matrix

# markov chain
markov <- function(x,P,n){
  seq <- x
  for(k in 1:(n-1)) {
    seq[k+1] <- sample(x, 1, replace=T, P[seq[k],])
  }
  return(seq)
}



hmmdat <- function(A, E, n) {
  observationset <- c(1:6)
  hiddenset <- c(1,2)
  x <- h <- matrix(NA, nr=n, nc=1)
  h[1] <- 1  # start with fair die
  x[1] <- sample(observationset,1,replace=T, E[h[1],])
  
  h <- markov(hiddenset, A, n)   # generate n hidden states
  
  # generate n observations 
  for(k in 1:(n-1)) {
    x[k+1] <- sample(observationset, 1, replace=T, E[h[k],])
  }
  
  out <- matrix(c(x,h), nrow=n, ncol=2, byrow=F)
  return(out)
}



E <- matrix(c(rep(1/6,6), rep(1/10,5), 1/2), 2, 6, byrow=T) # emission matrix
A <- matrix(c(0.95, 0.05, 0.1, 0.9),2,2,byrow=T)   # transition matrix

dat <- hmmdat(A,E,100)
colnames(dat) <- c("observation", "hidden state")
rownames(dat) <- 1:100
t(dat)
