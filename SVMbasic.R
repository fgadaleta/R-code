## SVM in R with kernlab package
library(kernlab)
library(ROCR)


#generate data

n <- 150  # number of data points
p <- 3    # dimensions

sigma <- 1 #variance of distrib.
meanpos <- 0   # mean of distrib of positive examples
meanneg <- 3   # mean of negative ex.
npos <- round(n/2)
nneg <- n - npos

#generate positive and negative examples
xpos <- matrix ( rnorm(npos*2, mean=meanpos, sd=sigma), npos, p)
xneg <- matrix ( rnorm(nneg*2, mean=meanneg, sd=sigma), nneg, p)
x <- rbind(xpos, xneg)

# generate labels
y <- matrix(c( rep(1,npos), rep(-1,nneg) ))


#visualize data
plot(x, col=ifelse(y>0,1,2))
legend("topleft", c('positive', 'negative'), col=seq(2), pch=1, text.col=seq(2))

#split training 80% and testing 20%
ntrain <- round(n*0.80)      # size of training set
tindex <- sample(n, ntrain)  #indeces of training samples
xtrain <- x[tindex,]         # training set 
xtest <- x[-tindex,]         # the rest is test set
ytrain <- y[tindex]
ytest <- y[-tindex]

istrain = rep(0,n)
istrain[tindex]=1
 
# visualize 
plot(x,col=ifelse(y>0,1,2), pch=ifelse(istrain==1,1,2))
legend("topleft",c('positive train','positive test','negative train','negative test'),col=c(1,1,2,2), pch=c(1,2,1,2), text.col=c(1,1,2,2))


# train SVM

svp <- ksvm(xtrain, ytrain, type="C-svc", kernel='rbf',C=100,scaled=c())
plot(svp, data=xtrain)

# predict labels on test
ypred = predict(svp, xtest)
table(ytest,ypred)  #contingency table
sum(ypred==ytest)/length(ytest) # accuracy

#compute at the prediction scores
ypredscore = predict(svp, xtest,type="decision")
table(ypredscore>0, ypred)

# ROC curves
pred <- prediction(ypredscore, ytest)

#plot ROC curve
perf <- performance(pred,measure="tpr",x.measure="fpr")
plot(perf)

# plot precision/recall curve
perf <- performance(pred, measure="prec", x.measure="rec")
plot(perf)

#plot accuracy as function of threshold
perf <- performance(pred, measure="acc")
plot(perf)