################################################
##  Manylasso Survey                          ##
##  (c) 2014 Francesco Gadaleta               ##    
##                                            ##  
################################################

library(plyr)  
source("importGold.R")
source("makegraph.R")
source("bench.R")
source("incremental.R")

library(amap)     # hcluster      
library(grpreg)   # group lasso 
library(hierNet)  # hierarchical lasso 
library(genlasso) # fused lasso 


## load data
data <- read.delim("../data/gnw_100_5/Ecoli-1_nonoise_multifactorial.tsv")
genes = names(data)  # keep right after data
p = ncol(data)
goldstd <- importGold('../data/gnw_100_5/Ecoli-1_goldstandard.tsv', labels=genes)



manylasso = function(data, scale=T, method=c("hier", "fused", "group", "lasso", "enet", "ridge", "onelasso", "oneenet", "oneridge")) {
  p<-ncol(data)
  n<-nrow(data)
  nets = list()
  
  timing = matrix(0, ncol=1, nrow = length(method))
  rownames(timing) = method
  colnames(timing) = "time[sec]"
  
  if(scale==T) { 
    data.scale <- scale(data)
    names(data.scale) = names(data)
  }
  
  for(i in 1:length(method)) {
    coeffs <- matrix(0, nrow = p, ncol = p)
    
    if(method[i] == "hier") {
      dbg <- sprintf("DBG %d, method %s lasso", i, method[i])
      print(dbg)
      
      start.time <- Sys.time()
      # for each covariate... 
      for(z in 1:p)  {
        j <- z  # advance progress bar
        #progress_bar_text$step()
      
        noti <- (1:p)[-z]
        yi <- data.scale[, z]
        xi <- data.scale[, noti]    
        
        fitcv = hierNet.path(xi,yi)
        singlefit = hierNet.cv(fitcv, xi,yi, nfolds=3)
        singlefit = hierNet(x= xi, y= yi, strong=T, lam=singlefit$lamhat, maxiter=1000)
        singlecoeffs = singlefit$bp - singlefit$bn # overall main effect estimated coefficients
        singleinter = round(singlefit$th,2)        # overall interactions estimated coefficients
        coeffs[z,-z] = singlecoeffs   # main effects are not averaged
        
        # additive interactions per gene
        coeffs[-z,-z] = coeffs[-z,-z] + ((singleinter)/p) # <pxp> matrix of interactions 
      }
      end.time <- Sys.time()
      timing["hier",1] <- difftime(end.time, start.time,units="secs")
      colnames(coeffs) = rownames(coeffs) = names(data)
      nets[[i]] = coeffs  
    }
   
    if(method[i] == "group") {
    groups = matrix(0, ncol = p, nrow = p)
    colnames(groups) = genes
    start.time <- Sys.time()
    
    # for each covariate... 
    for(z in 1:p)  {
      dbg <- sprintf("DBG regressing covariate %d", z)
      print(dbg)
      
      noti <- (1:p)[-z]
      yi <- data.scale[, z]
      xi <- data.scale[, noti]    
      
      # add correlated variables to the same group
      tmp = abs(cor(xi))
      t = hcluster(tmp)
      memb <- cutree(t, k = 3)  
      # TODO group differently (from goldstandard)
      #which(goldstd[1,] > 0, arr.ind = T)
      
      #fit = cv.grpreg(xi, yi, group=seq(1,p-1,by=1), nfolds=3)
      fit = cv.grpreg(xi, yi, group=memb, nfolds=10)
      cvlambda = fit$lambda.min
      fitgrp <- grpreg(xi,yi, lambda=cvlambda, group=memb, penalty = "grLasso")
      betas = as.matrix(fitgrp$beta)
      betas = betas[-1,]
      coeffs[z,-z] = betas 
    }
    end.time <- Sys.time()
    timing["group",1] <- difftime(end.time, start.time,units="secs") 
    colnames(coeffs) = rownames(coeffs) = names(data)
    nets[[i]] = coeffs  
    }
    
    if(method[i] == "fused"){
      start.time <- Sys.time()
      # for each covariate... 
      for(z in 1:p)  {
        dbg <- sprintf("DBG regressing covariate %d", z)
        print(dbg)
        
        noti <- (1:p)[-z]
        yi <- data.scale[, z]
        xi <- data.scale[, noti]    
        
        # add correlated variables to the same group
        tmp = abs(cor(xi))
        t = hcluster(tmp)
        numclusters = 10
        memb <- cutree(t, k = numclusters)  
        neigh = matrix(0, ncol = p-1, nrow = p-1)
        for(k in 1:numclusters) {
          fus = which(memb==k)
          neigh[fus, fus] = 1 
        }
        gr = graph.adjacency(neigh, mode = "undirected")
        
        fnet = fusedlasso(y=yi, X=xi, graph = gr, verbose=F)
        #fnet = fusedlasso1d(y=yi, X=xi, minlam=3, verbose=T)
        betas = as.matrix(fnet$beta)
        betas = betas[, dim(betas)[2]]
        coeffs[z,-z] = betas 
      }
      end.time <- Sys.time()
      timing["fused",1] <- difftime(end.time, start.time,units="secs")
      colnames(coeffs) = rownames(coeffs) = names(data)
      nets[[i]] = coeffs  
    }
    
    ## permutation-based methods
    if(method[i] == "lasso"){
      start.time <- Sys.time()
      netlasso = incrlasso(data=data, perm=p*5, repl=F, fanout=1, shrinkage=1, 
                       best=p, cutoff=1, cold=sample(1:p, p, replace=F), numcpu=4)
      end.time <- Sys.time()
      timing["lasso",1] <- difftime(end.time, start.time,units="secs")  
      nets[[i]] = netlasso
    }
    
    if(method[i] == "ridge"){
      start.time <- Sys.time()
      ridgelasso = incrlasso(data=data, perm=p*5, repl=F, fanout=1, shrinkage=0, 
                           best=p, cutoff=1, cold=sample(1:p, p, replace=F), numcpu=4)
      end.time <- Sys.time()
      timing["ridge",1] <- difftime(end.time, start.time,units="secs")
      nets[[i]] = ridgelasso
    }
    
    if(method[i] == "enet"){
      start.time <- Sys.time()
      enetlasso = incrlasso(data=data, perm=p*5, repl=F, fanout=1, shrinkage=0.5, 
                             best=p, cutoff=1, cold=sample(1:p, p, replace=F), numcpu=4)
      end.time <- Sys.time()
      timing["enet",1] <- difftime(end.time, start.time,units="secs")
      nets[[i]] = enetlasso
    }
    ## end of permutation-based 
 
    if(method[i] == "onelasso"){
      start.time <- Sys.time()
      # for each covariate... 
      for(z in 1:p)  {
        dbg <- sprintf("DBG regressing covariate %d", z)
        print(dbg)
        minlambda = 0
        
        noti <- (1:p)[-z]
        yi <- data.scale[, z]
        xi <- data.scale[, noti]    
        
        if(z == 1) {
          fit <- cv.glmnet(xi, yi)
          minlambda = fit$lambda.min
        }
        
        lnet = glmnet(x = xi, y = yi, lambda = minlambda, alpha = 1)
        betas = as.matrix(lnet$beta)
        coeffs[z,-z] = betas[, dim(betas)[2]] 
      }
      end.time <- Sys.time()
      timing["onelasso",1] <- difftime(end.time, start.time,units="secs")
      colnames(coeffs) = rownames(coeffs) = names(data)
      nets[[i]] = coeffs  
    }
    
    if(method[i] == "oneridge"){
      start.time <- Sys.time()
      # for each covariate... 
      for(z in 1:p)  {
        dbg <- sprintf("DBG regressing covariate %d", z)
        print(dbg)
        minlambda = 0
        
        noti <- (1:p)[-z]
        yi <- data.scale[, z]
        xi <- data.scale[, noti]    
        
        if(z == 1) {
          fit <- cv.glmnet(xi, yi)
          minlambda = fit$lambda.min
        }
        
        lnet = glmnet(x = xi, y = yi, lambda = minlambda, alpha = 0)
        betas = as.matrix(lnet$beta)
        coeffs[z,-z] = betas[, dim(betas)[2]] 
      }
      end.time <- Sys.time()
      timing["oneridge",1] <- difftime(end.time, start.time,units="secs")
      colnames(coeffs) = rownames(coeffs) = names(data)
      nets[[i]] = coeffs  
    }
    
    if(method[i] == "oneenet"){
      start.time <- Sys.time()
      # for each covariate... 
      for(z in 1:p)  {
        dbg <- sprintf("DBG regressing covariate %d", z)
        print(dbg)
        minlambda = 0
        
        noti <- (1:p)[-z]
        yi <- data.scale[, z]
        xi <- data.scale[, noti]    
        
        if(z == 1) {
          fit <- cv.glmnet(xi, yi)
          minlambda = fit$lambda.min
        }
        
        lnet = glmnet(x = xi, y = yi, lambda = minlambda, alpha = 0.5)
        betas = as.matrix(lnet$beta)
        coeffs[z,-z] = betas[, dim(betas)[2]] 
      }
      end.time <- Sys.time()
      timing["oneenet",1] <- difftime(end.time, start.time,units="secs")
      colnames(coeffs) = rownames(coeffs) = names(data)
      nets[[i]] = coeffs  
    }    
  }
  
  #save timing stats
  write.csv(timing, file="timing_manylasso.csv")
  
  names(nets) = method
  return(nets)
}


networks = manylasso(data)

olassonet = round(networks[["onelasso"]], 2)
oridgenet = round(networks[["oneridge"]], 2)
oenet     = round(networks[["oneenet"]], 2)

lassonet  = round(networks[["lasso"]], 2)
enet      = round(networks[["enet"]], 2)
ridgenet  = round(networks[["ridge"]], 2)

fusenet   = round(networks[["fused"]], 2)
groupnet  = round(networks[["group"]], 2)
hiernet   = round(networks[["hier"]], 2)


# zero out weak interactions at probs level
olassonet = abs(olassonet)
olassonet[olassonet < quantile(as.vector(olassonet), probs = 0.95)] = 0 

oridgenet = abs(oridgenet)
oridgenet[oridgenet < quantile(as.vector(oridgenet), probs = 0.95)] = 0 

oenet = abs(oenet)
oenet[oenet < quantile(as.vector(oenet), probs = 0.95) ] = 0 

groupnet = abs(groupnet)
groupnet[groupnet < quantile(as.vector(groupnet), probs = 0.95))] = 0 

lassonet = abs(lassonet)
lassonet[lassonet < quantile(as.vector(lassonet), probs = 0.95)] = 0 

fusenet = abs(fusenet)
fusenet[fusenet < quantile(as.vector(fusenet), probs = 0.95)] = 0 

ridgenet = abs(ridgenet)
ridgenet[ridgenet < quantile(as.vector(ridgenet), probs = 0.95)] = 0 

enet = abs(enet)
enet[enet < quantile(as.vector(enet), probs = 0.95)] = 0 

hiernet = abs(hiernet)
hiernet[hiernet < quantile(as.vector(hiernet), probs = 0.95)] = 0  

# compare with gold standard
res1 = bench(goldstd, fusenet) 
res2 = bench(goldstd, hiernet) 
res3 = bench(goldstd, groupnet) 

res4 = bench(goldstd, lassonet) 
res5 = bench(goldstd, ridgenet) 
res6 = bench(goldstd, enet) 

res7 = bench(goldstd, olassonet)
res8 = bench(goldstd, oridgenet)
res9 = bench(goldstd, oenet)

# collect and arrange results 
result = matrix(0, ncol = 10)
result = rbind(res1, res2, res3, res4, res5, res6, res7, res8, res9)
result = cbind(result, 0)
rownames(result) = c("fused", "hier", "group", "lasso", "ridge", "enet", "onelasso", "oneridge", "oneenet")
colnames(result) = c("gold", "pred", "tp","fp","tn","fn (missed)","mcc", "tpr", "fpr", "acc", "time")
timing = read.csv(file = "timing_manylasso.csv", header = T)
result[as.character(timing$X), ncol(result)] = timing$time.sec.

# write results
filename = sprintf("result_manylasso_%d_1.csv", dim(data)[2])
write.csv(result, file=filename)

