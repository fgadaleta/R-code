library(glmnet)
library(lars)

# Warning: if lasso, watch categorical variables

# n = 4100, p = 10000 (biallelic SNPS) 
# Format: [Individual ID; SNP1 allele1; SNP1 allele2; ..... SNP10,000 allele1; SNP10,000 allele2]
# read first nrows individuals
geno <- read.table("data/QTLMAS16/QTLMASXVI.gnt", header = FALSE, nrows = 3040)   

# Format: [individual; Trait1; Trait2; Trait3]
trait <- read.table("data/QTLMAS16/QTLMASXVI.dat", header = FALSE, nrows = 3000) 

#Format: [chromosome; SNP id; Position in cM; Position in pb;]
mapsnp <- read.table("data/QTLMAS16/QTLMASXVI.map", header = FALSE) 

# y <- response (trait pheno)
# x <- genotype of individual 
individuals <- trait[,1]   # individual indices

y <- trait[, 2]  

hist(trait[,4])

# set values of the first len SNPs (1 allele per time)
##idx1 <- seq(1, by = 2, len = 10000)
##idx1 <- idx1+1
##idx2 <- seq(0, by=2, len=10000) 
##idx2 <- idx2[-1]
##idx2 <- idx2 + 1


#x <- as.matrix(geno[individuals, idx1])
x <- as.matrix(geno[individuals, ])


# dominant model
# allele1 = 1, allele2 = 1   homozygote dominant   ==>  1
# allele1 = 1, allele2 = 2   hetero dominant       ==>  1
# allele1 = 2, allele2 = 1   hetero dominant       ==>  1
# allele1 = 2, allele2 = 2   homozygote recessive  ==> -1

numrows = dim(x)[1]
numsnps = 10000
snps = matrix(data=rep(1, numrows*numsnps), nrow=numrows, ncol=numsnps)


######################################################
# Dominant model
######################################################

for (i in 1:numrows){     # for each row
  for(j in 1:numsnps){    # for each col
    al1 <- 2*j
    al2 <- 2*j+1
    #snps[i,j] = 1   # default case is 1  take me out of the loop
    
    if ( (x[i, al1] == 2) & (x[i, al2] == 2) )   
      snps[i,j] = -1
    }
  }


######################################################
# TODO recessive model
######################################################


# alpha=1 is LASSO
# alpha=0 is RidgeRegression
fit <- glmnet(snps, y, alpha=0.1, family="gaussian")
plot(fit)
coefs <- coef(fit)
barplot(coefs[-1,1], main="QTL analysis", xlab="SNP position", ylab="Individual SNP effect")
deviance(fit)

#mae vs mse
cv.fit = cv.glmnet(snps,y, standardize=FALSE, alpha=0.1, type.measure="mse")
plot(cv.fit)
coefs <- coef(cv.fit, s="lambda.min")
#plot(coefs)
# coefs[-1,1]  takes coefficients of all variables but the intercept 
barplot(coefs[-1,1], main="QTL analysis", xlab="SNP position", ylab="Individual SNP effect")

train_pred <- predict(fit, x, type="coef")
compare = cbind(test, train_pred)

lfit <- lars(x=snps, y=y, trace=T, max.steps=1000, use.Gram=F)

# TODO

# map SNPs to genes

# build gene network from previous mapping

