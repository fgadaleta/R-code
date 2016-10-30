library(ggplot2)
qplot(data$pspF, data$rstB, data, geom=c("point", "smooth"))
qplot(data$pspF, data=data, fill=factor(data$rstB, labels=20))
qplot(data$pspF, data$rstB, data)

col = factor(data$pspF, levels=3)
qplot(data$pspF, data=data, geom="density", fill=col)
g <- ggplot(data, aes(data$pspF, data$rstB))
g + geom_point(color="steelblue", size=4, alpha=1/1.2) + theme_bw(base_family = "Times")
g + geom_point() + geom_smooth()

library(ggplot2)
g <- ggplot(movies, aes(votes, rating))
print(g)
g + geom_smooth() + geom_point()

airquality = transform(airquality, Month = factor(Month))
qplot(Wind, Ozone, data = airquality, facets = . ~ Month)

