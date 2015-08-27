####################################################################
## thr33some prj. expr-expr net construction                      ##
## Copyright 2014 Francesco Gadaleta                              ##
## University of Liege - Belgium                                  ##
####################################################################

require(snowfall)
library(MASS)

source("extract_utils.R")


#exprfile = "./data/GBM/GLIO_Gene_Expression.txt"
#methfile = "./data/GBM/GLIO_Methy_Expression.txt"

build_expr_methyl_net <- function(exprdata, methyldata, start=1, end=p, ncpu=4) 
{
  #methdata <- read.delim(methfile, header=T) 
  #exprdata <- read.delim(exprfile, header=T)
  namesexpr <- colnames(exprdata)       # samples expr
  namesmethyl <- colnames(methyldata)     # samples methyl
  p <- dim(exprdata)[2]
  
  assert(end <= p, "Not enough variables in your dataset!")
  
  # prepare adjacency matrix
  expr_methyl_net <- matrix(0, nrow = p, ncol = p)
  
  # prepare data for parallel computation 
  # expression - methylation data
  dbg <- sprintf("Preparing expr-methyl data for parallel computation")
  print(dbg)
  expr_meth_split = {}; j=1
  for(i in start:end){
    yi = as.numeric(as.matrix(exprdata[,i]))
    expr_meth_split[[j]] <- scale(as.matrix(cbind(yi, as.matrix(methyldata[,]))))
    j = j+1
    rm(yi)
  }
  
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
  
  sfInit(parallel = T, cpus = ncpu, type="SOCK")
  sfLibrary(glmnet)
  
  cvlambda = 0.1
  best = 4    # select the first best associations
  sfExport("cvlambda", "best")
  dbg <- sprintf("Inferring expr-methyl associations")
  print(dbg)
  fit_expr_meth <- sfLapply(expr_meth_split, singlefit)
  sfStop()
  rm(expr_meth_split)   # this is huge, free pls
  
  
  # name mapping
  singlemap_expr_methyl <- function(dat) {
    mnames = colnames(dat)
    selected = as.character(t(mnames[which(dat !=0, arr.ind=F)]))
    mappedgenenames = extract_gene_from_methyl(selected)    # selected genes from methyl positions
    selectedgenes = {}
    nmethyl = length(mappedgenenames)
    for(j in 1:nmethyl) {
      selectedgenes = c(selectedgenes, which(mappedgenenames[j] == namesexpr))
    }
    return(selectedgenes)
  }
  
  
  # prepare mapping data for parallel computation 
  mg_map = {}
  for(i in 1:length(fit_expr_meth)) {
    mg_map[[i]] = as.data.frame(t(fit_expr_meth[[i]]))
  }
  
  sfInit(parallel = T, cpus = 4, type="SOCK")
  sfExport("extract_gene_from_methyl", "namesexpr")
  dbg <- sprintf("Mapping methyl to gene space")
  print(dbg)
  res_expr_methyl <- sfLapply(mg_map, singlemap_expr_methyl)
  sfStop()
  
  # build adj matrix (expr-methyl)
  dbg <- sprintf("Building adjacency matrix expr-methyl")
  print(dbg)
  for(i in 1:length(res_expr_methyl)){
    selectedgenes = res_expr_methyl[[i]]
    expr_methyl_net[i, selectedgenes] = 1
  }
  
  # save expr-methyl network to disk 
  dbg <- sprintf("Saving expr-methyl network to disk")
  print(dbg)
  tmpfile = sprintf("./results/expr_meth/expr_methyl_subnet_%i_start_%i_end_%i.csv", p, start, end)
  expr_methyl_df = as.data.frame(expr_methyl_net)
  rownames(expr_methyl_df) = colnames(expr_methyl_df) = namesexpr
  write.table(x=expr_methyl_df, file=tmpfile)
  
  return(expr_methyl_df)
}