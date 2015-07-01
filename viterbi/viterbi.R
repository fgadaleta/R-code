# Viterbi algorithm
# 
# Run in Rstudio hmm_casino first, in order to get the observations in dat
# 

viterbi <- function(A,E,x) {
  v <- matrix(NA, nr=length(x), nc=dim(A)[1])
  v[1,] <- 0
  v[1,1] <- 1
  
  for(i in 2:length(x)) {
    for(l in 1:dim(A)[1]) {
      v[i,l] <- E[l,x[i]] * max(v[(i-1),]*A[l,])
        }
  }
  return(v)  
}

vit <- viterbi(A,E,dat[,1])
vitrowmax <- apply(vit, 1, function(x) which.max(x))  # select hidden state with max value
hiddenstate <- dat[,2]

table(hiddenstate, vitrowmax)  # performance of the algorithm 
datt <- cbind(dat, vitrowmax)
colnames(datt) <- c("observations", "hidden state", "predicted state")
t(datt);

sum(datt[,2]==datt[,3])
