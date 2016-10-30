library(gbm)
library(dplyr)

train = sample(1:699, 300, replace = F)
test = -train

data = BreastCancer
cs = factor(BreastCancer$Class)
classes = as.numeric(cs)

classes[classes == 2] = 0
data = cbind(data, classes)

?gbm

data.train = data[train,]
data.test = data[-train,]
model = gbm.fit(x=data.train[,2:10], y = data.train[,12], n.trees = 5000, distribution = "bernoulli", shrinkage = 0.001, verbose = T)
summary(model)
best.iter = gbm.perf(model)


f.predict <- predict(model, data.test[, c(2:10,12)], best.iter)
f.predict[f.predict>0] = 1
f.predict[f.predict<0] = 0

data.test[,12][1:10]
f.predict[1:10]
table(data.test[,12],f.predict)
