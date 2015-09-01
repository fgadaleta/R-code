##############################################
#  R Code for canonical correlation analysis #
##############################################

# We will use the built-in iris data set.
# We will consider the entire data set (all three species)
data(iris)
attach(iris)

# We will standardize the variables first
# by dividing by each column's standard deviation:
# (we will remove column 5, the species labels)

iris.std <- sweep(iris[,-5], 2, sqrt(apply(iris[,-5],2,var)), FUN="/")

sepal.meas <- iris.std[,1:2]
petal.meas <- iris.std[,3:4]

### Doing the CCA the long way:

# Finding blocks of the correlation matrix:

R11 <- cor(sepal.meas)
R22 <- cor(petal.meas)
R12 <- c(cor(sepal.meas[,1], petal.meas[,1]), cor(sepal.meas[,1], petal.meas[,2]),
         cor(sepal.meas[,2], petal.meas[,1]), cor(sepal.meas[,2], petal.meas[,2]))
R12 <- matrix(R12, ncol=ncol(R22), byrow=T) # R12 has q2 columns, same as number of petal measurements
R21 <- t(R12)  # R21=transpose of R12

# Finding the E1 and E2 matrices:

E1 <- solve(R11) %*% R12 %*% solve(R22) %*% R21
E2 <- solve(R22) %*% R21 %*% solve(R11) %*% R12

# print(E1)
# print(E2)

eigen(E1)
eigen(E2)

# The canonical correlations are:

canon.corr <- sqrt(eigen(E1)$values)
canon.corr

# The canonical variates are based on the eigenvectors of E1 and E2:

# a1 = (0.922, -0.388)
# b1 = (0.943, -0.333)
# a2 = (0.457, 0.890)
# b2 = (-0.679, 0.734)

# Only the first canonical correlation is really substantial:

# u1 = 0.92*Sepal.Length - 0.39*Sepal.Width
# v1 = 0.94*Petal.Length - 0.33*Petal.Width

# Plotting the first set of canonical variables:

u1 <- as.matrix(iris.std[,1:2]) %*% as.matrix(eigen(E1)$vectors[,1])
v1 <- as.matrix(iris.std[,3:4]) %*% as.matrix(eigen(E2)$vectors[,1])
plot(u1,v1)
cor(u1,v1)

# Plotting the second set of canonical variables:

u2 <- as.matrix(iris.std[,1:2]) %*% as.matrix(eigen(E1)$vectors[,2])
v2 <- as.matrix(iris.std[,3:4]) %*% as.matrix(eigen(E2)$vectors[,2])
plot(u2,v2)
cor(u2,v2)

### Doing CCA using the built-in cancor function:

cancor(sepal.meas, petal.meas)

# The canonical correlations are the same as the ones we found,
# The canonical variates are a little different because the cancor 
# function works with the centered data rather than the original data.