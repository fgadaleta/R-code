# SVM machines in R with e1071 package 
# (this package is  an interface to libsvm)

library('e1071')

data1 <- seq(1,10,by=2)
classes1 <- c('a','a','a','b','b')
test1 <- seq(1,10,by=2) + 1
model1 <- svm(data1,classes1,type='C', kernel='linear')
summary(model1)

pred <- predict(model1,data1)
#pred <- fitted(model1)
table(pred,classes1)


data2 <- seq(1,10)
classes2 <- c('b','b','b','a','a','a','a','b','b','b')
model2 <- svm(data2, classes2, type='C', kernel='radial')  
table(predict(model2,data2), classes2)   # with linear kernel not predicting right


bcdata <- read.csv(url("http://www.potschi.de/svmtut/breast-cancer-wisconsin.data"), head=TRUE)
names(bcdata)

#separate data from the classification
databcall <- subset(bcdata, select=c(-Samplecodenumber, -Class))
names(databcall)
classesbcall <- subset(bcdata, select=Class)

#training set
databctrain <- databcall[1:400,]
classesbctrain <- classesbcall[1:400,]


#test set
databctest <- databcall[401:699, ]
classesbctest <- classesbcall[401:699,]

# build SVM model
model <- svm(databctrain, classesbctrain)

# validate model (if not working here do not go ahead)
pred <- predict(model, databctest)
table(pred, t(classesbctest))

#improve model
obj <- tune(svm, train.x=databctrain, train.y=classesbctrain, validation.x=databctest,
     validation.y=classesbctest, ranges = list(gamma = 2^(-1:1),cost=2^(2:4)),
     control = tune.control(sampling = "fix"))



