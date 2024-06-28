
require(magic)

precissionMatrix3DFMMm <- function(paramArrayArg, t){
  d <- dim(paramArrayArg)[3]
  m <- dim(paramArrayArg)[1]
  
  alphas <- paramArrayArg[,3,1]
  omegas <- paramArrayArg[,5,1]
  deltas <- as.matrix(apply(paramArrayArg, 3, function(x){x[,2]*cos(x[,4])}))
  gammas <- as.matrix(apply(paramArrayArg, 3, function(x){-x[,2]*sin(x[,4])}))
  
  # Pre-calculus tl, u(alpha, omega, tl), (alpha, omega, tl)
  tl <- apply(cbind(alphas, omegas), 1, function(x){2*atan(x[2]*tan((t-x[1])/2))})
  u <- t(apply(tl, 1, sin))
  v <- t(apply(tl, 1, cos))
  if(nrow(u)<ncol(v)){u <- t(u); v <- t(v)}
  
  # Derivatives dMu/dParam
  blockM <- as.matrix(rep(1, length(t)*d))
  blockDelta <- u
  for (i in 1:(d-1)) {
    blockDelta <- adiag(blockDelta, u)
  }
  
  blockGamma <- v
  for (i in 1:(d-1)) {
    blockGamma <- adiag(blockGamma, v)
  }
  
  if(m==1){deltas <- t(deltas); gammas <- t(gammas)} # Cuidado las dimensiones
  
  blockAlphas <- sapply(1:m, function(x){(omegas[x] + (1 - omegas[x]^2)*(1-v[,x])/(2*omegas[x]))*(deltas[x, 1]*u[,x] - gammas[x, 1]*v[,x])})
  
  for (i in 1:(d-1)) {
    blockAlphas <- rbind(blockAlphas, sapply(1:m, function(x){(omegas[x] + (1 - omegas[x]^2)*(1-v[,x])/(2*omegas[x]))*(deltas[x, i]*u[,x] - gammas[x, i]*v[,x])}))
  }
  
  blockOmegas <- sapply(1:m, function(x){u[,x]*(-deltas[x, i]*u[,x] + gammas[x, i]*v[,x])/omegas[x]})
  for (i in 1:(d-1)) {
    blockOmegas <- rbind(blockOmegas, sapply(1:m, function(x){u[,x]*(-deltas[x, i]*u[,x] + gammas[x, i]*v[,x])/omegas[x]}))
  }
  
  F0 <- cbind(blockM, blockDelta, blockGamma, blockAlphas, blockOmegas)
  return(t(F0)%*%F0)
}

precissionMatrixFMM1 <- function(A, alpha, beta, omega, t){
  delta <- A*cos(beta)
  gamma <- -A*sin(beta)
  tl <- 2*atan(omega*tan((t-alpha)/2))
  u <- sin(tl)
  v <- cos(tl)
  derivM <- rep(1, length(t))
  derivDelta <- v
  derivGamma <- u
  derivAlpha <- (omega + (1 - omega^2)*(1-v)/(2*omega))*(delta*u - gamma*v)
  derivOmega <- u*(-delta*u + gamma*v)/omega
  F0 <- cbind(derivM, derivDelta, derivGamma, derivAlpha, derivOmega)
  return(t(F0)%*%F0)
}

precissionMatrixFMMm <- function(A, alpha, beta, omega, t){
  
  # Pre-calculus tl, u(alpha, omega, tl), (alpha, omega, tl)
  m <- length(alpha)
  delta <- A*cos(beta)
  gamma <- -A*sin(beta)
  tl <- apply(cbind(alpha, omega), 1, function(x){2*atan(x[2]*tan((t-x[1])/2))})
  u <- t(apply(tl, 1, sin))
  v <- t(apply(tl, 1, cos))
  if(nrow(u)<ncol(v)){u <- t(u); v <- t(v)}
  
  # Derivatives dMu/dParam
  derivM <- as.matrix(rep(1, length(t)))
  derivDelta <- v
  derivGamma <- u
  
  derivAlpha <- sapply(1:m, function(x){(omega[x] + (1 - omega[x]^2)*(1-v[,x])/(2*omega[x]))*(delta[x]*u[,x] - gamma[x]*v[,x])})
  derivOmega <- sapply(1:m, function(x){u[,x]*(-delta[x]*u[,x] + gamma[x]*v[,x])/omega[x]})
  
  F0 <- cbind(derivM, derivDelta, derivGamma, derivAlpha, derivOmega)
  return(t(F0)%*%F0)
}

