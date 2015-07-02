#############################################
#   (c) 2014 Francesco Gadaleta             #
#############################################


# ROC curves 
# comparison ROC - AUC

library(multtest)
library(ROCR)

data(golub)   # 3051 genes and 38 mRNA samples 
labels <- golub.cl

gol.true <- factor(golub.cl, levels=0:1, labels=c("ALL", "not ALL"))
gol.pred <- factor(golub[1042,] > 1.07, levels= c("TRUE", "FALSE"), labels=c("ALL", "not ALL"))

table(gol.pred, gol.true)


pred <- prediction(golub[1042,], gol.true)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, "auc")
acc <- performance(pred, "acc")
plot(acc)

plot(perf)

help(performance)
