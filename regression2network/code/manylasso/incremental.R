################################################
##  Manylasso Survey                          ##
##  (c) 2014 Francesco Gadaleta               ##    
##                                            ##  
################################################

library(glmnet)
library(snowfall)  
library(MASS)

## define p before calling
incrlasso <- function(data, scale=T, shrinkage=1, perm=100, 
                      cold=sample(1:p, p, replace=F), repl=T, fanout=1, 
                      best=2, cutoff, numcpu=4, start=1, end=p)
  
  {
  p<-ncol(data)
  n<-nrow(data)
  genenames <- names(data)
  numperm = perm
  coeffs <- matrix(0, nrow = p, ncol = p)
  #plambda = 0.001    # default lambda value (very small penalty)
  
  # randomly fill active list of nodes to be checked
  checked = rep(0, times=p)
  active = cold
  tbc = {}   # to be checked
  ccc = 0    # debug 
  discarded = 0 # number of covariates discarded by significance test
  
  if(scale==TRUE) { 
    data.scale <- scale(data)
    names(data.scale) = names(data)
  }
  
  ##compute lambda_CV
  lambda_cv = {} 
  numcov = {}
  lcv_max = round(p*.1)
  tst = sample(1:p, lcv_max, replace=F)
  print("Computing lambda_cv on 10% of covariates ...")
  for(t in 1:lcv_max) {
    noti <- (1:p)[-tst[t]]  
    ytest = data.scale[, tst[t]]     # response
    xtest = data.scale[, noti]       # covariates
    
    fit = cv.glmnet(x=xtest, y=ytest, standardize=F, intercept=F, alpha=shrinkage, type.measure="mae")  # mse    
    cvlambda = fit$lambda.min   # if cv.glmnet      
    nonzero = length(which(coef(fit)!= 0))
    
    lambda_cv = c(lambda_cv, cvlambda)
    numcov = c(numcov, nonzero)
    
    dbg <- sprintf("%d of %d", t,lcv_max)
    print(dbg)
  }  
  selected = median(numcov)
  plambda = mean(lambda_cv[which(numcov>=selected)])
  print("...done")
  
  ## set progress bar
  progress_bar_text <- create_progress_bar("text")
  progress_bar_text$init(p)  
  sfInit(parallel = TRUE, cpus = numcpu, type="SOCK") 
  
  # Function to be parallelised - takes loop index as argument.
  singleperm <- function(i)
  {
    #permutation of response
    res = sample(1:n, replace=repl)
    yinew = yi[res]
    #perturbation of response
    #vi = var(yi)
    #vi = 0.0001
    #yinew = yi + rnorm(n=1,mean=0, sd=vi)   #0.01
    #yinew = rnorm(n, 1, 0.01)
    fitnew = glmnet(x=xi, y=yinew, standardize=F, intercept=F, alpha=1, lambda=plambda) 
    singlecoeffs <- as.matrix(t(coef(fitnew, s="lambda.min")))
    singlecoeffs = singlecoeffs[,-1]                      # remove intercept
    ##fitnew = lars(xi, yinew)
    ##singlecoeffs <- predict(fitnew, s=rho, type="coefficients", mode="lambda")$coef #same rho is used for setting the weights and for lasso.
    singlecoeffs = which(singlecoeffs != 0 )
    singlecoeffs = names(singlecoeffs)
    pbetanew = {}
    
    # transform in index form    
    if(length(singlecoeffs) > 0 ) {
      for(t in 1:length(singlecoeffs)) {
        idx = match(singlecoeffs[t], names(data)) 
        pbetanew = c(pbetanew, idx)           # index of selected variable
      }  
    }
    return(pbetanew)
  }  
  
  # for each covariate... 
  # TODO we can use another dataset to select which covariate 
  # to discover interactions from how do we choose those? 
  # this can save a lot of computational effort
  
  while(length(active)>0)
  {    
    dbg <- sprintf("DBG %d", active)
    print(dbg)
    
    for(z in 1:length(active))
    {
      ccc = ccc+1
      j <- z  # advance progress bar
      progress_bar_text$step()
      
      noti <- (1:p)[-active[z]]
      yi <- data.scale[, active[z]]
      xi <- data.scale[, noti]    
      pbeta = {}     # vector of predicted beta 
      
      ## dbg
      dbg <- sprintf("regression of covariate %d (total checked %d)", active[z], ccc)
      print(dbg)
      
      fit = glmnet(xi, yi, standardize=F, intercept=F, alpha=shrinkage, lambda=plambda) 
      singlecoeffs <- as.matrix(t(coef(fit, s="lambda.min")))
      singlecoeffs = singlecoeffs[,-1]                      # remove intercept
      
      maincoeffsidx = sort(abs(singlecoeffs), decreasing=T, index.return=T)$ix
      effects = singlecoeffs[maincoeffsidx]
      maineffects = which(effects != 0 )
      
      if(best <= length(maineffects))
        maineffects = names(maineffects)[1:best] 
      else
        maineffects = names(maineffects)
      
      singlecoeffs = maineffects
      
      # transform in index form and save to pbeta (predicted beta)
      if(length(singlecoeffs) > 0 ) {
        for(t in 1:length(singlecoeffs)) {
          idx = match(singlecoeffs[t], names(data))
          pbeta = c(pbeta, idx)
        }  
      }
      
      # add all neighbors, to be checked next
      checked[active[z]] = 1    # this covariates has been checked
      ##tbc = c(tbc, pbeta)       # these covariates have to be checked next 
      
      #############################################
      #permutation-based significance test
      #############################################
      # counter of selections from permutated data
      counterpar = rep(0, times=p)
      # Make data available to other R instances / nodes
      sfExport(list = c("xi", "yi", "n", "plambda"))
      # To load a library on each R instance / node
      sfClusterEval(library(glmnet))
      # Use a parallel RNG to avoid correlated random numbers
      sfClusterSetupRNG()
      coefs <- sfClusterApplyLB(1:numperm, singleperm) 
      
      ## merge all slaves into counter and everything stays the same from here on
      # merge results from single slave
      for(t in 1:length(coefs)) {
        c <- unlist(coefs[t])
        counterpar[c] = counterpar[c] + 1
      }
      
      frequency = counterpar/numperm   # TODO this is P(connection(a,b))
      counter = counterpar
      
      # return selected coefficients from significance testing
      rem = {}
      #nsel = fanout  # should vary according to the p-value
      
      if(length(pbeta)>0) {
        sigcoefs <- pbeta
        sigcount <- counter[sigcoefs] 
        # increasingly order the counter (we want the least selected after permutation)
        coefidx <- sort(sigcount, decreasing=F, index.return=T)$ix  
        
        print("preselected coefficients")
        print(pbeta)
        
        print("counters of preselected coefficients")
        print(sigcount)  
        
        f = factor(sigcount)
        selfan = {}
        
        # add the first <fanout> factors
        for(fan in 1:fanout) {
          sel = which(sigcount==levels(f)[fan], arr.ind=T)
          ###selfan = c(pbeta[sel], selfan)
          themin = min(sigcount)
          
          # empirical distribution of counters 
          if(length(sigcount)>2) {
            ed <- density(sigcount)
            pr_min <- sum(ed$y[ed$x<themin])/(sum(ed$y))
          }
          else {
            pr_min = 0
          }
          
          dbg <- sprintf("prob(select %d for %d times) = %f", pbeta[sel], min(sigcount), pr_min)
          print (dbg)
          
          if(pr_min < cutoff) {
            selfan = c(pbeta[sel], selfan)    
          }
          else {discarded = discarded + 1}
        }
        
        sigcoefs = selfan
       
        dbg <- sprintf("+ %d ", sigcoefs)
        print (dbg)
        tbc = c(tbc, sigcoefs)       # these covariates have to be checked next 
       
        coeffs[active[z], sigcoefs] = counter[sigcoefs]             # adds to matrix of coefs 
      }
      
      if(length(pbeta) == 0) {
        coeffs[active[z],] = 0  
      }    
      
      ### save tmp results to file
      tmpfile = sprintf("tmp_subnet_nodes_%i_perm_%i_%i.csv", dim(data)[2], perm, ccc)
      write.matrix(format(coeffs, scientific=FALSE), 
                   file = paste("../results/", tmpfile, sep="/"), sep=",")
      
    }   # for all vars in active
    
    
    tbc = tbc[which(checked[tbc] == 0)] 
    tbc = unique(tbc)        # remove duplicates
    
    print("Checking next")
    dbg <- sprintf("%d ", tbc)
    print(dbg)
    
    active = tbc
  }   # end of while 
  
  sfStop()  # we can stop the cluster now
  cat("\n")   # progress bar
  
  dbg <- sprintf("%d covariates discarded by significance test", discarded)
  print (dbg)
  
  colnames(coeffs) = rownames(coeffs) = genenames
  return(coeffs)
}

