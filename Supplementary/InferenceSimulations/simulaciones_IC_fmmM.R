library(openxlsx)
library(R.utils)
sourceDirectory("funcionesFMM/")

between <- function(x, inf, sup){x>inf && x<sup}

precissionMatrixFMM <- function(A, alpha, beta, omega, t){
  
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

logL <- function(params, m, vData, t, varErr){
  M <- params[1]
  delta <- params[1+(1:m)]
  gamma <- params[1+(1:m)+m]
  alpha <- params[1+(1:m)+2*m]
  omega <- params[1+(1:m)+3*m]
  
  yFit <- M + apply(apply(rbind(delta, gamma, alpha, omega), 2, function(x, t){
        x[1]*cos(2*atan(x[4]*tan((t-x[3])/2))) + x[2]*sin(2*atan(x[4]*tan((t-x[3])/2)))}, t = t), 1, sum)
  
  return(-sum(0.5*(1/varErr)*(vData-yFit)^2))
}

m <- 5
n <- 200

M <- 416
A <- c(63.46, 121.9, 420.18, 117.22, 127.93)
alpha <- c(4.9, 5.61, 5.75, 5.89, 1.3)
beta <- c(3.37, 1.24, 3.33, 5.39, 3.14)
omega <- c(0.13, 0.03, 0.03, 0.03, 0.19)

delta <- A*cos(beta)
gamma <- -A*sin(beta)
t <- seqTimes(n)
nreps <- 20
maxiter <- 40
sigma <- 35
IN1 <- matrix(0, nrow = nreps, ncol = 4*m)
IN1B <- matrix(0, nrow = nreps, ncol = 4*m)
IN2 <- matrix(0, nrow = nreps, ncol = 4*m)
IN3 <- matrix(0, nrow = nreps, ncol = 4*m)

estimadores <- matrix(0, nrow = nreps, ncol = 4*m)
SEF0Sub <- matrix(0, nrow = nreps, ncol = 4*m)
SEF0SubB <- matrix(0, nrow = nreps, ncol = 4*m)
SEF0 <- matrix(0, nrow = nreps, ncol = 4*m)
SEHess <- matrix(0, nrow = nreps, ncol = 4*m)

for (i in 1:nreps) {
  ff <- generateFMM(M = M, A = A, alpha = alpha, beta = beta, omega = omega, 
                    sigmaNoise = sigma, length.out = n, plot = F)
  
  # ML ESTIMATION
  exECG <- fitFMM(ff$y, nback = 5, parallelize = T, maxiter = maxiter)
  plotFMM(exECG)
  plotFMM(exECG, components = T)
  sum(exECG@R2)
  sigmaHat <- sqrt(exECG@SSE/(n-(1+4*m)))
  waveOrder <- order((exECG@alpha+pi)%%(2*pi))
  Ahat <- exECG@A[waveOrder]
  alphahat <- exECG@alpha[waveOrder]
  betahat <- exECG@beta[waveOrder]
  omegahat <- exECG@omega[waveOrder]
  optFMM <- optim(c(exECG@M, Ahat*cos(betahat), -Ahat*sin(betahat), 
                    alphahat, omegahat), 
                  fn = logL, control = list(fnscale = -1), t = t, m = m, 
                  vData = ff$y, varErr = sigmaHat^2, hessian = T)
  
  trueParameters <- c(delta, gamma, alpha, omega)
  

  out <- try({
    sdsB <- precissionMatrixFMM(Ahat[2:4], alphahat[2:4], betahat[2:4], omegahat[2:4], t)
    sdsB <- sqrt(diag(sigmaHat^2*solve(sdsB)))[-1]
  }, silent = T)
  
  if(!class(out) == 'try-error'){
    estimatedParameters <- optFMM$par[-1]
    estimadores[i,] <- estimatedParameters
    
    ############################# F0 en submodelos ###############################
    
    sds <- apply(X = rbind(Ahat, alphahat, betahat, omegahat), MARGIN = 2, function(x, t){
      varMat <- sigmaHat^2*solve(precissionMatrixFMM(x[1], x[2], x[3], x[4], t))
      return(sqrt(diag(varMat)))}, t = t)
    sds <- sds[-1,]
    SE <- as.vector(t(sds))
    
    lInf <- estimatedParameters - 2*SE
    lSup <- estimatedParameters + 2*SE
    namesICs <- paste0( rep(c("delta", "gamma","alpha", "omega"), each = 5), 
                        rep(c("P", "Q", "R", "S", "T"), 4))
    
    ICs <- rbind(lInf, trueParameters, lSup)
    colnames(ICs) <- namesICs
    paramIn <- apply(ICs, 2, function(x){between(x[2], x[1], x[3])})
    
    IN1[i,] <- paramIn
    SEF0Sub[i,] <- SE
    
    ############################ F0 en submodelos ################################
    ## 1-P 2-QRS 3-T
    sdsB <- precissionMatrixFMM(Ahat[2:4], alphahat[2:4], betahat[2:4], omegahat[2:4], t)
    sdsB <- sqrt(diag(sigmaHat^2*solve(sdsB)))[-1]
    SEB <- SE
    SEB[2:4] <- sdsB[1:3]
    SEB[7:9] <- sdsB[4:6]
    SEB[12:14] <- sdsB[7:9]
    SEB[17:19] <- sdsB[10:12]
    
    lInfB <- estimatedParameters - 2*SEB
    lSupB <- estimatedParameters + 2*SEB
    
    ICsB <- rbind(lInfB, trueParameters, lSupB)
    colnames(ICsB) <- namesICs
    paramInB <- apply(ICsB, 2, function(x){between(x[2], x[1], x[3])})
    
    IN1B[i,] <- paramInB
    SEF0SubB[i,] <- SEB
    
    ############################# F0 COMPLETA ####################################
    
    # Con F0 COMPLETA
    tF0F0 <- precissionMatrixFMM(Ahat, alphahat, betahat, omegahat, t = t)
    sds2 <- sigmaHat^2*solve(tF0F0)
    SE2 <- sqrt(diag(sds2))[-1]
    
    lInf2 <- estimatedParameters - 2*SE2
    lSup2 <- estimatedParameters + 2*SE2
    
    ICs2 <- rbind(lInf2, trueParameters, lSup2)
    colnames(ICs2) <- namesICs
    paramIn2 <- apply(ICs2, 2, function(x){between(x[2], x[1], x[3])})
    
    IN2[i,] <- paramIn2
    SEF0[i,] <- SE2
    ################################## OPTIM #####################################
    
    SE3 <- sqrt(diag(solve(-optFMM$hessian)))[-1]
    
    lInf <- estimatedParameters - 2*SE3
    lSup <- estimatedParameters + 2*SE3
    
    ICs3 <- rbind(lInf, trueParameters, lSup)
    colnames(ICs3) <- namesICs
    paramIn3 <- apply(ICs3, 2, function(x){between(x[2], x[1], x[3])})
    
    IN3[i,] <- paramIn3
    SEHess[i,] <- SE3
  }else{
    print("!!")
    i = i-1
  }
  
}
namesICs <- paste0( rep(c("delta", "gamma","alpha", "omega"), each = 5), 
                    rep(c("P", "Q", "R", "S", "T"), 4))
# Coverturas

Covert_F0Sub <- apply(IN1, 2, sum, na.rm = T)
Covert_F0SubB <- apply(IN1B, 2, sum, na.rm = T)
Covert_F0 <- apply(IN2, 2, sum, na.rm = T)
Covert_Hess <- apply(IN3, 2, sum, na.rm = T)
Covertures <- (rbind(Covert_F0Sub,Covert_F0SubB, Covert_F0, Covert_Hess)*100)/nreps
colnames(Covertures) <- namesICs

# Sesgo y varianza de los estimadores

sesgo <- round(apply(apply(estimadores, 1, function(x){x-trueParameters}), 1, mean), 4)
names(sesgo) <- namesICs
varianzas <- apply(estimadores, 2, var)
sd <- round(sqrt(varianzas), 5)

names(sd) <- namesICs
sesgoSdEstimadores <- rbind(sesgo, sd)

# Media y sd de los errores estándar

meanSE_F0Sub <- round(apply(SEF0Sub, 2, mean, na.rm = T), 5)
sdSE_F0Sub <- round(sqrt(apply(SEF0Sub, 2, var, na.rm = T)), 5)
meanSE_F0SubB <- round(apply(SEF0SubB, 2, mean, na.rm = T), 5)
sdSE_F0SubB <- round(sqrt(apply(SEF0SubB, 2, var, na.rm = T)), 5)
meanSE_F0 <- round(apply(SEF0, 2, mean, na.rm = T), 5)
sdSE_F0 <- round(sqrt(apply(SEF0, 2, var, na.rm = T)), 5)
meanSE_Hess <- round(apply(SEHess, 2, mean, na.rm = T), 5)
sdSE_Hess <- round(sqrt(apply(SEHess, 2, var, na.rm = T)), 5)

hatSE <- rbind(meanSE_F0Sub, sdSE_F0Sub, meanSE_F0SubB, sdSE_F0SubB, meanSE_F0, sdSE_F0, meanSE_Hess, sdSE_Hess)
colnames(hatSE) <- namesICs

# Salidas

Covertures
sesgoSdEstimadores
hatSE

# Guardamos todos 

resumen <- as.data.frame(rbind(Covertures, sesgoSdEstimadores, hatSE))
resumen
# write.xlsx(resumen, "AnalisisInfFMM5-500.xlsx",colNames = TRUE, rowNames = TRUE, 
#            colWidths = 10)

# ests <- as.data.frame(estimadores)
# colnames(ests) <- namesICs
# 
# write.xlsx(ests, "estimadoresFMM5-500.xlsx",colNames = TRUE, rowNames = TRUE, 
#            colWidths = 10)

# colnames(SEF0Sub) <- namesICs
# colnames(SEF0) <- namesICs
# 
# write.xlsx(as.data.frame(cbind(SEF0Sub, SEF0)), "ErroresEstandar-500.xlsx", colNames = T, 
#            colWidths = 7)

d=1
resumenDeltas <- apply(resumen[1:4,1:(m*d)], 1, mean)
resumenGammas <- apply(resumen[1:4,(m*d+1):(2*m*d)], 1, mean)
resumenAlphas <- apply(resumen[1:4,(m*2*d+1):(m*2*d+m)], 1, mean)
resumenOmegas <- apply(resumen[1:4,(m*2*d+m+1):(m*2*d+2*m)], 1, mean)
cobertMedias <- as.data.frame(cbind(resumenDeltas, resumenGammas, resumenAlphas, resumenOmegas))
# write.xlsx(cobertMedias, "AnalisisInfFMM5-500cobMedias.xlsx", colNames = TRUE, rowNames = TRUE,
#            colWidths = 10)
d = 1
resumenP <- apply(resumen[1:4,5*(1:(2*d+2))-4], 1, mean)
resumenQ <- apply(resumen[1:4,5*(1:(2*d+2))-3], 1, mean)
resumenR <- apply(resumen[1:4,5*(1:(2*d+2))-2], 1, mean)
resumenS <- apply(resumen[1:4,5*(1:(2*d+2))-1], 1, mean)
resumenT <- apply(resumen[1:4,5*(1:(2*d+2))], 1, mean)
resumenWaves  <- as.data.frame(cbind(resumenP, resumenQ, resumenR, resumenS, resumenT))
# write.xlsx(resumenWaves, "AnalisisInfFMM5-500cobWaves.xlsx", colNames = TRUE, rowNames = TRUE,
#            colWidths = 10)
















