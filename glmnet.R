###########################################################
# GLMNET
###########################################################
library(glmnet)

#read data
train <- read.csv("data/train.csv", header=TRUE)
test <- read.csv("data/test.csv", header=TRUE)
subset <- 1000

trainGLMNET <- train[1:subset,]
trainGLMNET$Target_Leaderboard = NULL
targetGLMNET <- trainset$Target_Leaderboard

#clean up the test set for prediction
testGLMNET <- test[1:subset,]
testGLMNET$Target_Leaderboard = NULL
testGLMNET$case_id = NULL
testGLMNET$train = NULL
testGLMNET$Target_Evaluate = NULL
testGLMNET$Target_Practice = NULL

#get the value of lambda through cross validation
mylambda <- cv.glmnet(as.matrix(trainGLMNET),
                      trainGLMNET[,1],family="multinomial",
                      type="auc",nfolds=10)
plot(mylambda,ylim=c(0,1))

best.lambda  <- mylambda$lambda.min

#build the model using that value of lambda
glmnet_model <- glmnet(as.matrix(trainGLMNET),
                       trainGLMNET[,1],family="multinomial",
                       lambda=best.lambda)

#predict
train_GLMNET <- predict(glmnet_model,type="response",as.matrix(trainGLMNET))
test_GLMNET <- predict(glmnet_model,type="response",as.matrix(testGLMNET))


########################################
#Generate a file for submission
########################################
testID  <- testset$case_id
predictions <- test_GLMNET
submit_file = cbind(testID,predictions)
write.csv(submit_file, file="GLM_benchmark.csv", row.names = FALSE)