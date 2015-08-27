####################################################################
### the "thr33some" project 
### (C) 2014 Francesco Gadaleta University of Liege 
####################################################################

require(snowfall)
library(MASS)
####################################################################
###                   load GBM raw data
GLIO_Methy_Expression <- read.delim("./data/GBM/GLIO_Methy_Expression.txt", header=F)
GLIO_Gene_Expression <- read.delim("./data/GBM/GLIO_Gene_Expression.txt", header=F)
GLIO_Survival <- read.delim("./data/GBM/GLIO_Survival.txt")
###                   load GBM expression
expr = as.data.frame(t(GLIO_Gene_Expression))

# reduce data dimension for testing only!!
expr = expr[, 1:100]

indexpr = as.character(expr[-1,1])  # samples 
namesexpr = {}                      # names of genes from expression dataset
# some tricky code to get the name of each gene 
namesline = expr[1,-1]
exprdata = expr[-1,-1]              # raw expression data 
exprdata= exprdata[-dim(exprdata)[1], ]   # remove NAs 
p = length(namesline)
for(i in 1:p)
{
  nameidx = length(levels(namesline[,i]))
  namesexpr = c(namesexpr, levels(namesline[[i]])[nameidx]) 
}
rm(namesline, expr)
### end of tricky code

###                   load GBM methylation
methyl = as.data.frame(t(GLIO_Methy_Expression))
indmethyl = as.character(methyl[-1,1])
namesmethyl = methyl[1,-1]
methyldata = methyl[-1,-1]          # raw methyl data
methyldata = methyldata[-dim(methyldata)[1],]
rm(methyl)
###                   load survival 
survdata = GLIO_Survival
survdatacases = survdata[which(survdata$Death==1),]
survcasesindex = which(survdata$Death==1, arr.ind=T)
survcontrolindex = which(survdata$Death==0, arr.ind=T)
####################################################################


####################################################################
# data extraction utilities
####################################################################

# one-to-one mapping methyl id into gene space
# usage: out = extract_gene_from_methyl(namesmethyl)
extract_gene_from_methyl <- function(namesmethyl)
{
  end = length(namesmethyl)
  genespace = {}
  for(i in 1:end) {
    ## get methyl_id
    methyl_id = namesmethyl[i]
    methylgene = strsplit(methyl_id, "_")[[1]][1]
    genespace = c(genespace, methylgene)
  }
  return(genespace)
}

gm = extract_gene_from_methyl(namesmethyl)
####################################################################

####################################################################
## data analysis utilities
####################################################################

## (for each gene regress methylation against expression)
## connect current gene to methylation-mapped genes and build gene-gene network
library(glmnet)

# prepare adjacency matrices
expr_methyl_net <- matrix(0, nrow = p, ncol = p)
expr_expr_net   <- matrix(0, nrow = p, ncol = p)


# prepare data for parallel computation 
# expression - methylation data
dbg <- sprintf("Preparing expr-methyl data for parallel computation")
print(dbg)
expr_meth_split = {}
for(i in 1:p){
  yi = as.numeric(as.matrix(exprdata[,i]))
  #xi = as.matrix(methyldata[,])  
  expr_meth_split[[i]] <- as.matrix(cbind(yi, as.matrix(methyldata[,])))
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

sfInit(parallel = T, cpus = 4, type="SOCK")
sfLibrary(glmnet)
# compute cv.lambda
#fit = cv.glmnet(y=as.numeric(as.matrix(exprdata[,1])), 
#                x=as.matrix(methyldata), 
#                alpha=1, intercept=F)
#cvlambda = fit$lambda.min
# TODO cross validation lambda 
cvlambda = 0.1
best = p    # select the first best associations
sfExport("cvlambda", "best")
dbg <- sprintf("Inferring associations expr-methyl")
print(dbg)
fit_expr_meth <- sfLapply(expr_meth_split, singlefit)
sfStop()
rm(expr_meth_split) # this is huge

#name mapping
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
#names = as.data.frame(namesmethyl)
for(i in 1:p) {
  mg_map[[i]] = as.data.frame(t(fit_expr_meth[[i]]))
  #names = as.data.frame(namesmethyl)
  #smap = rbind(tmp, namesmethyl)
  #mg_map[[i]] = smap
  #rm(smap)
}
#rm(names)

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
tmpfile = sprintf("./results/expr_methyl_net_%i.csv", p)
expr_methyl_df = as.data.frame(expr_methyl_net)
rownames(expr_methyl_df) = colnames(expr_methyl_df) = namesexpr
write.table(x=expr_methyl_df, file=tmpfile)



############ expr-expr net #############
# scale data
exprdata.scale <- scale(as.matrix(data.matrix(exprdata[survcasesindex,])))


# prepare data for parallel computation 
# expression - methylation data
dbg <- sprintf("Preparing expr-expr data for parallel computation")
print(dbg)
expr_expr_split = {}
for(i in 1:p){
  yi = as.numeric(as.matrix(exprdata.scale[,i]))
  xi = as.matrix(exprdata.scale[,-i])  
  expr_expr_split[[i]] <- as.matrix(cbind(yi, xi))
  rm(xi, yi)
}

sfInit(parallel = T, cpus = 4, type="SOCK")
sfLibrary(glmnet)

# compute cv.lambda
# TODO check some and average
cvlambda = 0
for(i in 1:round(p*0.1)){
  ncv = sample(x=1:p, 1, replace=F)
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
rm(expr_expr_split)   # this is huge 


# prepare mapping data for parallel computation 
dbg <- sprintf("Preparing data for expr-expr mapping")
print(dbg)
mg_map = {}
for(i in 1:p) {
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


sfInit(parallel = T, cpus = 4, type="SOCK")
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
tmpfile = sprintf("./results/expr_expr_net_%i.csv", p)
expr_expr_df = as.data.frame(expr_expr_net)
rownames(expr_expr_df) = colnames(expr_expr_df) = namesexpr
write.table(x=expr_expr_df, file=tmpfile)


## visualisation 
library(igraph)
g <- graph.adjacency(expr_expr_net, weighted=T, mode="undirected", diag=F)
V(g)$name = namesexpr   # fix this when generalising
png(filename="./results/expr_expr_net.png", width=1000, height=1000)
plot.igraph(g, vertex.size = 5)
dev.off()

g <- graph.adjacency(expr_methyl_net, weighted=T, mode="undirected", diag=F)
V(g)$name = namesexpr   # fix this when generalising
png(filename="./results/expr_methyl_net.png", width=1000, height=1000)
plot.igraph(g, vertex.size = 5)
dev.off()


## TODO 
## expr <- SNP 
## (for each gene regress SNPs against expression)  
## connect current gene to SNP-mapped genes and build gene-gene network
