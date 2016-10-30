# makes the random forest submission

library(randomForest)

train <- read.csv("data/train.csv", header=TRUE)
test <- read.csv("data/test.csv", header=TRUE)

labels <- as.factor(train[,1])
train <- train[,-1]

rf <- randomForest(train, labels, xtest=test, ntree=100)
predictions <- levels(labels)[rf$test$predicted]

plot(rf)

write(predictions, file="rf_benchmark.csv", ncolumns=1) 


names(train)


library(party)
train <- read.csv("data/train.csv", header=TRUE)
test <- read.csv("data/test.csv", header=TRUE)

#ctree and cforest expect labels to be factors
train$label <- as.factor(train$label)
class(train$label)

partree <- ctree(label~., data=train[1:1000,])
plot(partree)

data.controls <- cforest_unbiased(ntree=100, mtry=28)

subset = 5000
cf<- cforest(label~., data=train[1:subset,], 
             control = data.controls)

table(predict(cf), train[1:subset,1])
sum(predict(cf) == train[1:subset,1])/subset


tr <- treeresponse(cf, newdata = test[1:100,])
