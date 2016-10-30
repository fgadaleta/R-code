library(e1071)
library(multtest)

data(golub)
gol.fac <- factor(golub.cl, levels=0:1, labels=c("ALL", "AML"))

svmest <- svm(t(golub), gol.fac, type= "C-classification", kernel="linear")
svmpred <- predict(svmest, t(golub))
table(svmpred, gol.fac)

# or this way
pred <- fitted(svmest)
table(pred, gol.fac)
