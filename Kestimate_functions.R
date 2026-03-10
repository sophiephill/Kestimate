
#####################################################################
####### Estimate K

## Solve Ax = b where A = U s V'
## if reg_params > 0 includes L2 regularization
tikh <- function(U, s, V, b, reg_params = 0) {
  utb <- t(U) %*% b
  zeta <- drop(s * utb)
  
  ll <- length(reg_params)
  xlam <- matrix(0, nrow = length(s), ncol = ll)
  
  for (i in 1:ll) {
    xlam[,i] <- drop(V %*% (zeta / (s^2 + reg_params[i])) )
    
  }
  colnames(xlam) <- reg_params
  return(list(coef = xlam) )
}

rmse <- function(x) sqrt(mean(x * x))
norm <- function(x) sqrt(sum(x * x))

force0 <- function(v) ifelse(v < 0, 0, v)
