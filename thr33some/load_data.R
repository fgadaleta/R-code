####################################################################
### the "thr33some" project 
### (C) 2014 Francesco Gadaleta University of Liege 
####################################################################
# load data
# scale data
# save data 

library(MASS)

####################################################################
###                   load GBM raw data
GLIO_Gene_Expression <- read.delim("./data/GBM/GLIO_Gene_Expression.txt", header=T)
GLIO_Methy_Expression <- read.delim("./data/GBM/GLIO_Methy_Expression.txt", header=T)
GLIO_Survival <- read.delim("./data/GBM/GLIO_Survival.txt", header=T)

###                   load GBM expression
exprdata = as.data.frame(t(GLIO_Gene_Expression))
# reduce data dimension for testing purposes only!!
exprdata = exprdata[, 1:500]
indexpr = rownames(exprdata[-1,])        # samples
namesexpr = colnames(exprdata)           # gene names
exprdata= exprdata[-dim(exprdata)[1], ]  # remove NAs 
p = length(namesexpr)                    # number of genes
#########################################################

###                   load GBM methylation
methyldata = as.data.frame(t(GLIO_Methy_Expression))
indmethyl = rownames(methyldata)              # samples
namesmethyl = colnames(methyldata)            # methyl names
methyldata = methyldata[-dim(methyldata)[1],] # remove NAs
#########################################################


###                   load survival 
survdata = GLIO_Survival
survdatacases = survdata[which(survdata$Death==1),]
survcasesindex = which(survdata$Death==1, arr.ind=T)
survcontrolindex = which(survdata$Death==0, arr.ind=T)




############ expr-expr net #############
# scale data
#exprdata.scale <- scale(as.matrix(data.matrix(exprdata[survcasesindex,])))

### save scaled data to disk
#tmpfile = sprintf("expression_data_scaled_samples_%i_genes_%i.csv", 
#                  dim(exprdata.scale)[1], dim(exprdata.scale)[2])
#write.matrix(format(exprdata.scale, scientific=FALSE), 
#             file = paste("./data/", tmpfile, sep="/"), sep=",")

####################################################################
