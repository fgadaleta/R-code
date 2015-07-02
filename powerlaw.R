#############################################
#   (c) 2014 Francesco Gadaleta             #
#############################################

library(gsl) 
library(numDeriv)
library(ggplot2)

data = sort(sim, decreasing=T)
qplot(data, geom="density", size=2, fill=factor(data, levels=seq(from=1,to=40,by=5))
xmins = unique(data)
dat = numeric(length(xmins)) 
z = sort(data)

for (i in 1:length(xmins)){
  xmin = xmins[i] # choose next xmin candidate
  z1 = z[z>=xmin] # truncate data below this xmin value
  n = length(z1)
  a = 1+ n*(sum(log(z1/xmin)))^-1 # estimate alpha using direct MLE
  cx = (n:1)/n
  cf = (z1/xmin)^(-a+1) 
  dat[i] = max(abs(cf-cx)) 
}

D = min(dat[dat>0],na.rm=TRUE) # find smallest D value xmin = xmins[which(dat==D)] # find corresponding xmin value
z = data[data>=xmin]
z = sort(z)
n = length(z)
alpha = 1 + n*(sum(log(z/xmin)))^-1 # get corresponding alpha estimate




# the following code, up to the rpowerlaw, 
#came from this website: http://www.rickwash.com/papers/cscw08-appendix/powerlaw.R
dpowerlaw <- function(x, alpha=2, xmin=1, log=F) 
  { 
  if (log)
  log(alpha-1) - log(xmin) - alpha * log(x / xmin) 
  else
    ((alpha - 1) / xmin) * ((x / xmin) ^ (-alpha)) 
}


ppowerlaw <- function(q, alpha=2, xmin=1, lower.tail=T, log.p = F) 
  { 
  p <- (q / xmin) ^ (- alpha + 1)                                                                     
  if (lower.tail)                                                                  
    p <- 1-p 
  if (log.p)
    p <- log(p) 
  p
}


qpowerlaw <- function(p, alpha=2, xmin=1, lower.tail=T, log.p = F) 
  { 
  if (!lower.tail)
  p <- 1-p if (log.p)
    p <- exp(p)
  xmin * ((1 - p) ^ (-1 / (alpha - 1)))
}


rpowerlaw <- function(n, alpha=2, xmin=1) { 
  qpowerlaw(runif(n, 0, 1), alpha, xmin)
}

testresult = numeric(2500) 

for (i in 1:2500){
  power = rpowerlaw(length(z),alpha,xmin) # parameters we found
  #randomly generate power law data using the
  w = ks.test(z,power) #using KS test to see how good the fit is if (w$p.value > 0.10){
  testresult[i] = 1
}

if (w$p.value <= 0.10){
  testresult[i] = 0} 
}

sum(testresult)



# this does the same automagically
library(poweRlaw)
m_pl = displ$new(floor(sim))
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl)
bs = bootstrap(m_pl, no_of_sims=5000, threads=10)
