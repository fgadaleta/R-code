library(ElemStatLearn)

data(galaxy)
?galaxy

# spectral decomposition
v=var(galaxy)

e = eigen(v, symmetric=T)
e$vectors %*% t(e$vectors)

e$vectors %*% diag(e$values) %*% t(e$vectors)
G = e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
G %*% G

Ginv = e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)
G %*% Ginv


W = sweep(galaxy, 2, colMeans(galaxy))
Z = as.matrix(W) %*% Ginv

# now original data have been standardised (mean=0, var=1)
colMeans(Z)
var(Z)


# compute Mahalanobis distance for each observation
# on standardised data
d = sqrt(apply(Z*Z, 1, sum))

#or let R implement all the above :)
f = sqrt(mahalanobis(galaxy, colMeans(galaxy), var(galaxy) ))

hist(d)
galaxy[d>4, ]   # observations at a distance greater than 4

plot(galaxy$radial.position, galaxy$velocity, pch=19, col=4)
points(galaxy$radial.position[d>4], galaxy$velocity[d>4], pch=19, col=2)





