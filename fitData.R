#################################################
# Automatic distribution fitting and selection
#
# Copyright 2014 worldofpiggy.com
#################################################
library(MASS)

#Usage: 
#data, numeric vector of observations of unknown distribution
#fit, a list of distributions to fit data
#sample, rate of subsampling (0.5 means that a sample 50% of data will be considered) 
fitData <- function(data, fit="gamma", sample=0.5){
  distrib = list()
  numfit <- length(fit)
  results = matrix(0, ncol=5, nrow=numfit)
  
  for(i in 1:numfit){
    if((fit[i] == "gamma") | 
         (fit[i] == "poisson") | 
         (fit[i] == "weibull") | 
         (fit[i] == "exponential") |
         (fit[i] == "logistic") |
         (fit[i] == "normal") | 
         (fit[i] == "geometric")
    ) 
      distrib[[i]] = fit[i]
    else stop("Provide a valid distribution to fit data" )
  }
  
  # take a sample of dataset
  n = round(length(data)*sample)
  data = sample(data, size=n, replace=F)
  
  for(i in 1:numfit) {
    if(distrib[[i]] == "gamma") {
      gf_shape = "gamma"
      fd_g <- fitdistr(data, "gamma")
      est_shape = fd_g$estimate[[1]]
      est_rate = fd_g$estimate[[2]]
   
      ks = ks.test(data, "pgamma", shape=est_shape, rate=est_rate)
      
      # add to results
      results[i,] = c(gf_shape, est_shape, est_rate, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "poisson"){
      gf_shape = "poisson"
      fd_p <- fitdistr(data, "poisson")
      est_lambda = fd_p$estimate[[1]]
      
      ks = ks.test(data, "ppois", lambda=est_lambda)
      # add to results
      results[i,] = c(gf_shape, est_lambda, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "weibull"){
      gf_shape = "weibull"
      fd_w <- fitdistr(data,densfun=dweibull,start=list(scale=1,shape=2))
      est_shape = fd_w$estimate[[1]]
      est_scale = fd_w$estimate[[2]]
      
      ks = ks.test(data, "pweibull", shape=est_shape, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_shape, est_scale, ks$statistic, ks$p.value)   
    }
    
    else if(distrib[[i]] == "normal"){
      gf_shape = "normal"
      fd_n <- fitdistr(data, "normal")
      est_mean = fd_n$estimate[[1]]
      est_sd = fd_n$estimate[[2]]
      
      ks = ks.test(data, "pnorm", mean=est_mean, sd=est_sd)
      # add to results
      results[i,] = c(gf_shape, est_mean, est_sd, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "exponential"){
      gf_shape = "exponential"
      fd_e <- fitdistr(data, "exponential")
      est_rate = fd_e$estimate[[1]]
      ks = ks.test(data, "pexp", rate=est_rate)
      # add to results
      results[i,] = c(gf_shape, est_rate, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "logistic"){
      gf_shape = "logistic"
      fd_l <- fitdistr(data, "logistic")
      est_location = fd_l$estimate[[1]]
      est_scale = fd_l$estimate[[2]]
      ks = ks.test(data, "plogis", location=est_location, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_location, est_scale, ks$statistic, ks$p.value)     
    }
  }
  results = rbind(c("distribution", "param1", "param2", "ks stat", "ks pvalue"), results)
  #print(results)
  return(results)
}


# let's try this out

#example: normal (easily confused by an 'equivalent' gamma)
data = rnorm(100000, mean=5, sd=0.75)
res = fitData(data,  fit=c("logistic", "normal", "exponential", "poisson"), sample=1)

#example: exponential
data = rexp(100000,rate=2.3)
res = fitData(data,  fit=c("logistic", "normal", "exponential", "poisson"), sample=1)


# example: gamma
data= rgamma(100000, shape=3.2, rate=1.3)
res = fitData(data,  fit=c("logistic", "normal", "exponential", "gamma", "poisson"), sample=0.1)
