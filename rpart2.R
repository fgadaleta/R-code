# gene selection from golub data

library(multtest)
library(rpart)

data(golub)
row.names(golub) <- paste("gene", 1:3051, sep="")
goldata <- data.frame(t(golub[1:3051,]))          # this is how data should be (covariates as column)
gol.fac <- factor(golub.cl, levels=0:1, labels=c("ALL", "AML"))
gol.rp <- rpart(gol.fac ~ ., data=goldata, method="class", cp=0.001)
plot(gol.rp, branch=0, margin=0.1); text(gol.rp, digits=3, use.n=TRUE) # this shows gene896
golub.gnames[896,]   # this is the best ALL-AML predictor 
