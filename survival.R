library(survival)
library(KMsurv)

data(aids)
attach(aids)
colnames(aids)
infect
induct
detach(aids)

data(psych)
attach(psych)
my.surv.obj <- Surv(age, age+time, death)
fit <- survfit(my.surv.obj~sex, data = psych)
detach(psych)

summary(fit)$surv
summary(fit)$time
summary(fit)$n.risk
summary(fit)$n.event

plot(fit, main="KM estimate", xlab = "time", ylab = "survival funct", 
     col = c("red", "blue"))



data(tongue)
my.surv <- Surv(tongue$time, tongue$delta)
fitt <- survfit(my.surv~tongue$type, conf.int = 0.95)
plot(fitt)
