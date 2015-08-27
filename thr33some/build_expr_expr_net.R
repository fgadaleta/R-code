####################################################################
## thr33some prj. expr-expr net construction                      ##
## Copyright 2014 Francesco Gadaleta                              ##
## University of Liege - Belgium                                  ##
####################################################################

require(snowfall)
library(MASS)

build_expr_expr_net <- function(exprdata, start=1, end=p, ncpu=4) 
  {
  data <- scale(as.matrix(data.matrix(exprdata)))
  data = data[-1,]    # remove X line (data dependent)
  namesexpr <- colnames(data)
  p <- dim(data)[2]
  
  # prepare adjacency matrix
  expr_expr_net   <- matrix(0, nrow = p, ncol = p)
  
  # single fitting function for parallel computation
  singlefit <- function(dat) {
    class(dat) = "numeric"
    fit = glmnet(y=dat[,1], x=dat[,-1], alpha=1, intercept=F, lambda=cvlambda)
    singlecoeffs <- as.matrix(t(coef(fit)))
    singlecoeffs <- singlecoeffs[,-1]  # remove intercept
    
    ## keeps only the <best> best 
    numnonzero = length(which(singlecoeffs!=0))
    if(best < numnonzero) { 
      maineffects = sort(abs(singlecoeffs), decreasing=T, index.return=T)$ix
      maineffects = maineffects[1:best]
      # sets the rest to 0
      singlecoeffs[-maineffects] = 0 
    }
    
    return(singlecoeffs) 
  }
  
  
  # prepare data for parallel computation 
  # expression-expression data
  dbg <- sprintf("Preparing expr-expr data for parallel computation")
  print(dbg)
  expr_expr_split = {}; j=1
  for(i in start:end){
    yi = as.numeric(as.matrix(data[,i]))
    xi = as.matrix(data[,-i])  
    expr_expr_split[[j]] <- as.matrix(cbind(yi, xi))
    j = j+1
    rm(xi, yi)
  }
  
  sfInit(parallel = T, cpus = ncpu, type="SOCK")
  sfLibrary(glmnet)
  
  # compute cv.lambda
  # TODO check some and average
  cvlambda = 0
  for(i in 1:round(p*0.1)){
    ncv = sample(x=1:length(expr_expr_split), 1, replace=F)
    fit = cv.glmnet(y=expr_expr_split[[ncv]][,1], 
                    x=expr_expr_split[[ncv]][,-1], 
                    alpha=1, intercept=F)
    cvlambda = cvlambda + fit$lambda.min
  }
  cvlambda = cvlambda/ncv
  cvlambda = 0.1
  best = 4    # select the first best associations
  sfExport("cvlambda", "best")
  dbg <- sprintf("Inferring expr-expr associations")
  print(dbg)
  fit_expr_expr <- sfLapply(expr_expr_split, singlefit)
  sfStop()
  rm(expr_expr_split)   # this is huge, free pls
  
  # prepare mapping data for parallel computation 
  dbg <- sprintf("Preparing data for expr-expr mapping")
  print(dbg)
  mg_map = {}
  
  for(i in 1:length(fit_expr_expr)) {
    tmp = as.data.frame(t(fit_expr_expr[[i]]))   # selected genes from fitting
    smap = rbind(tmp, namesexpr[-i])
    mg_map[[i]] = smap
    rm(smap)
  }
  
  # maps selected genes onto genespace 
  # return: index of selected gene 
  singlemap_expr_expr <- function(dat) {
    mnames = dat[2,]
    selected = t(mnames[which(dat[1,] !=0, arr.ind=T)])
    selectedgenes = {}
    nexpr = length(selected)
    for(j in 1:nexpr) {
      selectedgenes = c(selectedgenes, which(selected[j] == namesexpr))
    }
    return(selectedgenes)
  }
  
  
  sfInit(parallel = T, cpus = ncpu, type="SOCK")
  sfExport("namesexpr")
  dbg <- sprintf("Mapping expr-expr selections")
  print(dbg)
  res_expr_expr <- sfLapply(mg_map, singlemap_expr_expr)
  sfStop()
  
  # build adj matrix (expr-expr)
  dbg <- sprintf("Building adjacency expr-expr matrix")
  print(dbg)
  for(i in 1:length(res_expr_expr)){
    selectedgenes = res_expr_expr[[i]]
    expr_expr_net[i, selectedgenes] = 1
  }
  
  # save expr-expr network to disk 
  dbg <- sprintf("Saving expr-expr network to disk")
  print(dbg)
  #FIXME create directory if not exist
  tmpfile = sprintf("./results/expr_expr/expr_expr_subnet_%i_start_%i_end_%i.csv", 
                    p, start,end)
  expr_expr_df = as.data.frame(expr_expr_net)
  rownames(expr_expr_df) = colnames(expr_expr_df) = namesexpr
  write.matrix(format(expr_expr_df, scientific=F), 
               file = tmpfile, sep=",")
  
  return(expr_expr_df)
}
