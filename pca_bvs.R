# A comparison between Principal Component Analysis and 
# Bayesian Variable Selection on the Golub (1999) dataset


library(multtest)
data(golub)         # expression of 3051 genes from 38 patients
y <- golub.cl       # classes of patients (ALL=0 AML=1)
gol.fac <- factor(class, levels=0:1, labels= c("ALL", "AML"))  # factors


golub[1042, gol.fac=="ALL"]
golub[1042, gol.fac=="AML"]
grep("CD33",golub.gnames[,2])   #index of gene CD33 in the golub dataset


sum(eigen(cor(golub))$values[1:5])/38*100 # coverage of total amount of variance in data
-eigen(cor(golub))$vec[,1:2][,1]          # weigths of first 2 eigenvectors for all genes


pca <- princomp(golub, center=TRUE, cor=TRUE, scores=TRUE)
o <- order(pca$scores[,2])
golub.gnames[o[1:10], 2]      # names of first 10 representative genes wrt 2nd component
golub.gnames[o[3041:3051],2]  # names of the last 10 genes wrt to 2nd component

#######################################################
# variable selection K&M
# Issue: we are forcing the model to select about 20% of variables 
# (which could lead to a Lindley's paradox situation)
library(rjags)

p <- dim(golub)[1]
n <- dim(golub)[2]

data=list(y=y,x=t(golub),n=n,p=p)
init=list(tau=1,alpha=0,betaT=rep(0,p),ind=rep(0,p))
modelstring="
  model {
    for (i in 1:n) {
      mean[i]<-alpha+inprod(x[i,],beta)
      y[i]~dnorm(mean[i],tau)
    }
    for (j in 1:p) {
      ind[j]~dbern(0.2)
      betaT[j]~dnorm(0,0.001)
      beta[j]<-ind[j]*betaT[j]
    }
    alpha~dnorm(0,0.0001)
    tau~dgamma(1,0.001)
  }
"
model=jags.model(textConnection(modelstring),
                 data=data,inits=init)
update(model,n.iter=500)
output=coda.samples(model=model,
                    variable.names=c("alpha","beta","ind","tau"),
                    n.iter=1000,thin=1)
print(summary(output))