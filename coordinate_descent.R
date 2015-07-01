# t(x)*x = t(x)*y
solve(t(x)%*%x, t(x)%*%y)
p = 2
i = 1 
i = 2

for(k in 1:10){
 for(i in 1:2) 
  theta[i] = t(x[,i])%*%(y - x[,-i]%*%t(theta[-i]))/ (t(x[,i])%*%x[,i])
}


theta[1] = 0
theta[2] = 0
