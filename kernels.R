#############################################
#   (c) 2014 Francesco Gadaleta             #
#############################################

library(kernlab)

data(spam)
dt <- as.matrix(spam[c(10:20, 3000:3010), -58])

rbf <- rbfdot(sigma = 0.5)

#rbf(o1,o2)
#exp(0.5 * (2 * crossprod(o1, o2) - crossprod(o1) - crossprod(o2)))
#exp(0.5 * (2 * t(o1)%*%o2 - t(o1)%*%o1 - t(o2)%*%o2))

K = kernelMatrix(kernel = rbf, x = dt)
yt <- as.matrix(as.integer(spam[c(10:20, 3000:3010), 58]))
yt[yt==2] = -1

# compute kernel expansion
kernelMult(rbf, dt, , yt)



# 1. example of ksvm
data(promotergene)
tindex <- sample(1:dim(promotergene)[1], 5)
genetrain <- promotergene[-tindex, ]
genetest <- promotergene[tindex, ]

gene <- ksvm(Class ~ ., data = genetrain, kernel = "rbfdot",
            kpar = "automatic", C = 60, cross = 3, prob.model = TRUE)
predict(gene, genetest, type="probabilities")
tmp <- genetest[1, -1]
#tmp["V2"] = 'a'
predict(gene, tmp)

# 2. example of ksvm
x <- rbind(matrix(rnorm(120), , 2), matrix(rnorm(120, mean = 3), , 2))
y <- matrix(c(rep(1, 60), rep(-1, 60)))
vp <- ksvm(x, y, type = "C-svc")
plot(vp, data=x)
newdata = matrix(c(2,0.1, 3,3), ncol = 2)
predict(vp, newdata)

# relevance vector machine
rvmm <- rvm(x,y, kernel="rbfdot", kpar=list(sigma=0.1))
yhat <- predict(rvmm, x)
plot(yhat)


# ranking algorithm
data(spirals)
ran <- spirals[rowSums(abs(spirals) < 0.55) == 2, ]
ranked <- ranking(ran, 54, kernel = "rbfdot", 
                  kpar = list(sigma = 100), edgegraph = TRUE)
ranked[54,2] <- max(ranked[-54,2])
c <- 1:86
op <- par(mfrow = c(1, 2), pty = "s")
plot(ran)
plot(ran, cex = c[ranked[, 3]]/40)


# online learning
x <- rbind(matrix(rnorm(90), , 2), matrix(rnorm(90) + 3, , 2))
y <- matrix(c(rep(1, 45), rep(-1, 45)), , 1)
on <- inlearn(2, kernel = "rbfdot", kpar = list(sigma = 0.2), type = "classification")
ind <- sample(1:90, 90)

for (i in ind) 
  on <- onlearn(on, x[i, ], y[i], nu = 0.03, lambda = 0.1)
sign(predict(on, x))

spc <- specc(x, centers=4)
plot(spc@centers)
plot(x)




