### the arbitrary constant
const <- 5

### function with two equations
twoFuncs <- function(par){
  x.1 = par^2
  x.2 = 2*(const - par)^2
  return(sum(c(x.1, x.2)))
}

### function with three equations
threeFuncs <- function(par){
  x.1 = par[1]^2
  x.2 = 2*par[2]^2
  x.3 = 3*(const - sum(par))^2
  return(sum(c(x.1, x.2, x.3)))
}


### function with four equations
fourFuncs <- function(par){
  x.1 = par[1]^2
  x.2 = 2*par[2]^2
  x.3 = 3*par[3]^2
  x.4 = -5*(const - sum(par))
  return(sum(c(x.1, x.2, x.3, x.4)))
}

### try constraining a different parameter
### should get same result as above
fourFuncs2 <- function(par){
  x.1 = par[1]^2
  x.2 = 2*(const - sum(par))^2
  x.3 = 3*par[2]^2
  x.4 = -5*par[3]
  return(sum(c(x.1, x.2, x.3, x.4)))
}


solution <- optim(0, twoFuncs, method="BFGS")
solution <- optim(numeric(2), threeFuncs, method="BFGS")
solution <- optim(numeric(3), fourFuncs, method="BFGS")
solution <- optim(numeric(3), fourFuncs2, method="BFGS")


solution$par



