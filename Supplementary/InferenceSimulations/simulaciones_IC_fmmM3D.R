
library(openxlsx)
library(dplyr)
library(R.utils)
library(magic)
sourceDirectory("funcionesFMM/")
source("auxMultiFMMtest2.R")

normBeat <- read.xlsx("param_config.xlsx", "NORM", )

leadParams <- normBeat[,-1]
n = 200
m = 5
d = 8
leadNames <- normBeat$Lead[1:d]
compNames <- (c("P", "Q", "R", "S", "T"))[1:m]
nreps = 200 
maxiter = 50
commonSigma <- 35

# Cargar desde aqui hasta la 104

{
  between <- function(x, inf, sup){x>inf && x<sup}
  precissionMatrixFMM3d <- function(paramArrayArg, t){
    d <- dim(paramArrayArg)[3]
    m <- dim(paramArrayArg)[1]
    
    alphas <- paramArrayArg[,3,2]
    omegas <- paramArrayArg[,5,2]
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
  namesDeltas <- paste0("delta_",rep(compNames, d), "_",rep(leadNames, each = m))
  namesGammas <- paste0("gamma_",rep(compNames, d), "_",rep(leadNames, each = m))
  namesAlphas <- paste0("alpha_",compNames)
  namesomega <- paste0("omega_",compNames)
  fullNames <- c(namesDeltas, namesGammas, namesAlphas, namesomega)
}

paramArray <- array(0, c(m, 5, d))
estimatesArray <- array(0, c(m, 5, d))

if(d==3){
  selectedLeads <- c(1,2,4)
}else if(d==5){
  selectedLeads <- c(1,2,4,6,7)
}else{
  selectedLeads <- 1:d
}

for (k in 1:d){
  paramArray[,,k] <- matrix(as.numeric(leadParams[selectedLeads[k],1:(5*m)]), ncol = 5, byrow = T)
}

# paramArray[4,3,] <- 5.91
alphas <- paramArray[,3,1]
omegas <- paramArray[,5,1]
As <- paramArray[,2,]
betas <- paramArray[,4,]
deltas <- c(As*cos(betas))
gammas <- c(-As*sin(betas))
trueParameters <- c(deltas, gammas, alphas, omegas)
t <- seqTimes(n)

################################################################################
# Estructuras
numCols <- (1+d)*2*m
IN1 <- matrix(0, nrow = nreps, ncol = numCols)
IN1B <- matrix(0, nrow = nreps, ncol = numCols)
IN2 <- matrix(0, nrow = nreps, ncol = numCols)
IN3 <- matrix(0, nrow = nreps, ncol = numCols)

estimadores <- matrix(0, nrow = nreps, ncol = numCols)
SEF0Sub <- matrix(0, nrow = nreps, ncol = numCols)
SEF0SubB <- matrix(0, nrow = nreps, ncol = numCols)
SEF0 <- matrix(0, nrow = nreps, ncol = numCols)
SEHess <- matrix(0, nrow = nreps, ncol = numCols)

################################################################################

#SSETrace <- 1:maxiter
#sesgosEjemplo <- matrix(0, ncol = 10, nrow = maxiter)

for (rep in 1:nreps) {

  cat("---------------------",rep,"----------------------\n\n")
  
  # Simulation
  simulatedData <- apply(paramArray, 3, function(x){generateFMM(x[1,1], x[,2], x[,3], x[,4], x[,5], 
                                                                sigmaNoise = commonSigma, length.out = n, plot = F)$y})
  colnames(simulatedData) <- leadNames
  
  # Fitting an estimated values
  paramsPerSignal <- fitMultiFMM(vDataMatrix = simulatedData, maxIter = maxiter, nBack = m,
                                 omegaMin = 0.001, omegaMax = 0.5, commonOmega = TRUE, weightError = FALSE)
  
  alphasHat <- paramsPerSignal[[1]][,3]
  waveOrder <- order((alphasHat+3*pi/4)%%(2*pi))
  for (k in 1:d){
    paramsPerSignal[[k]] <- paramsPerSignal[[k]][waveOrder,]
    estimatesArray[,,k] <- as.matrix(paramsPerSignal[[k]][,-6])
  }
  
  alphasHat <- alphasHat[waveOrder]; omegasHat <- estimatesArray[,5,1]
  AsHat <- estimatesArray[,2,]; betasHat <- estimatesArray[,4,]
  deltasHat <- c(As*cos(betasHat)); gammasHat <- c(-As*sin(betasHat))
  estimatedParameters <- c(deltasHat, gammasHat, alphasHat, omegasHat)
  names(estimatedParameters) <- fullNames
  estimadores[rep,] <- estimatedParameters
  
  r2 <- sapply(1:length(paramsPerSignal), function(x) sum(paramsPerSignal[[x]][substr(rownames(paramsPerSignal[[x]]),1,1)!="X",]$Var))
  sigmaHat <- sqrt(sum((1-r2)*apply(simulatedData, 2, function(vData){sum((vData - mean(vData))^2)}))/(d*n - (1+(1+d)*2*m)))
  # sigmaHat2 <- sqrt(sum((1-r2)*apply(simulatedData, 2, function(vData){sum((vData - mean(vData))^2)}))/(d*n))
  
  ############################# F0 en submodelos ###############################
  
  sds <- apply(estimatesArray, MARGIN = 1, function(x, t){
    x2 <- array(x, c(1,5,d))
    varMat <- sigmaHat^2*solve(precissionMatrixFMM3d(x2, t))
    return(sqrt(diag(varMat)))}, t = t)
  
  sds <- sds[-1,]
  SE <- as.vector(t(sds))
  
  lInf <- estimatedParameters - 2*SE
  lSup <- estimatedParameters + 2*SE
  
  ICs <- rbind(lInf, trueParameters, lSup)
  colnames(ICs) <- fullNames
  paramIn <- apply(ICs, 2, function(x){between(x[2], x[1], x[3])})
  IN1[rep,] <- paramIn
  SEF0Sub[rep,] <- SE
  
  ############################ F0 en submodelos ################################
  ## 1-P 2-QRS 3-T
  
  sdsB <- precissionMatrixFMM3d(estimatesArray[c(2:4),,], t)
  sdsB <- sqrt(diag(sigmaHat^2*solve(sdsB)))[-1]
  SEB <- SE
  SEB[rep(0:(d-1), each = 3)*m + 2:4] <- sdsB[1:(3*d)]
  SEB[rep(0:(d-1), each = 3)*m + 2:4 + d*5] <- sdsB[(3*d+1):(2*3*d)]
  SEB[(2*m*d+2):(2*m*d+4)] <- sdsB[(2*3*d+1):(2*3*d+3)]
  SEB[(2*m*d+7):(2*m*d+9)] <- sdsB[(2*3*d+4):(2*3*d+6)]

  lInfB <- estimatedParameters - 2*SEB
  lSupB <- estimatedParameters + 2*SEB

  ICsB <- rbind(lInfB, trueParameters, lSupB)
  colnames(ICsB) <- fullNames
  paramInB <- apply(ICsB, 2, function(x){between(x[2], x[1], x[3])})

  IN1B[rep,] <- paramInB
  SEF0SubB[rep,] <- SEB

  ############################# F0 COMPLETA ####################################
  
  tF0F0 <- precissionMatrixFMM3d(estimatesArray, t)
  seHat2 <- sqrt(diag(sigmaHat^2*solve(tF0F0)))
  seHat2 <- seHat2[-1]
  
  lInf2 <- estimatedParameters - 2*seHat2
  lSup2 <- estimatedParameters + 2*seHat2
  
  ICs2 <- rbind(lInf2, trueParameters, lSup2)
  colnames(ICs2) <- fullNames
  paramIn2 <- apply(ICs2, 2, function(x){between(x[2], x[1], x[3])})
  IN2[rep,] <- paramIn2
  SEF0[rep,] <- seHat2

}

Covert_F0Sub <- apply(IN1, 2, sum, na.rm = T)
Covert_F0SubB <- apply(IN1B, 2, sum, na.rm = T)
Covert_F0 <- apply(IN2, 2, sum, na.rm = T)
Coverages <- (rbind(Covert_F0Sub,Covert_F0SubB, Covert_F0)*100)/nreps
colnames(Coverages) <- fullNames

##### END

sesgo <- apply(apply(estimadores, 1, function(x){x-trueParameters}), 1, mean)
names(sesgo) <- fullNames
seHat <- apply(estimadores, 2, sd)
aux <- rbind(trueParameters, sesgo, seHat)
sesgo


meanlIC_F0Sub <- 4*apply(SEF0Sub, 2, mean)
sdSE_F0Sub <- sqrt(apply(SEF0Sub, 2, var))
meanlIC_F0SubB <- 4*apply(SEF0SubB, 2, mean)
sdSE_F0SubB <- sqrt(apply(SEF0SubB, 2, var))
meanlIC_F0 <- 4*apply(SEF0, 2, mean)
sdSE_F0 <- sqrt(apply(SEF0, 2, var))
hatSE <- rbind(meanlIC_F0Sub, sdSE_F0Sub, meanlIC_F0SubB, sdSE_F0SubB, meanlIC_F0, sdSE_F0)
colnames(hatSE) <- fullNames
hatSE


resumen <- as.data.frame(rbind(Coverages, sesgo, hatSE))
resumenDeltas <- apply(resumen[1:3,1:(m*d)], 1, mean)
resumenGammas <- apply(resumen[1:3,(m*d+1):(2*m*d)], 1, mean)
resumenAlphas <- apply(resumen[1:3,(m*2*d+1):(m*2*d+m)], 1, mean)
resumenOmegas <- apply(resumen[1:3,(m*2*d+m+1):(m*2*d+2*m)], 1, mean)
cobertMedias <- as.data.frame(cbind(resumenDeltas, resumenGammas, resumenAlphas, resumenOmegas))

cobertMedias


resumen[,(m*2*d+1):(m*2*d+2*m)] # Solo saca alphas y omegas
cobertMedias

# write.xlsx(resumen, "AnalisisInfFMM3D-M3-500.xlsx", colNames = TRUE, rowNames = TRUE,
#           colWidths = 10)
# write.xlsx(cobertMedias, "AnalisisInfFMM3D-M3-500cobMedias.xlsx", colNames = TRUE, rowNames = TRUE,
#            colWidths = 10)
# estimadores <- data.frame(estimadores)
# colnames(estimadores) <- fullNames
# write.xlsx(estimadores, "AnalisisInfFMM3D-M3-500estimadores.xlsx",colNames = TRUE, rowNames = TRUE,
#            colWidths = 10)


################################################################################

resumenDeltas <- apply(resumen[1:3,1:(m*d)], 1, mean)
resumenGammas <- apply(resumen[1:3,(m*d+1):(2*m*d)], 1, mean)
resumenAlphas <- apply(resumen[1:3,(m*2*d+1):(m*2*d+m)], 1, mean)
resumenOmegas <- apply(resumen[1:3,(m*2*d+m+1):(m*2*d+2*m)], 1, mean)
cobertMedias <- as.data.frame(cbind(resumenDeltas, resumenGammas, resumenAlphas, 
                                    resumenOmegas))
cobertMedias
################################################################################

resumenP <- apply(resumen[1:3,5*(1:(2*d+2))-4], 1, mean)
resumenQ <- apply(resumen[1:3,5*(1:(2*d+2))-3], 1, mean)
resumenR <- apply(resumen[1:3,5*(1:(2*d+2))-2], 1, mean)
resumenS <- apply(resumen[1:3,5*(1:(2*d+2))-1], 1, mean)
resumenT <- apply(resumen[1:3,5*(1:(2*d+2))], 1, mean)

resumenWaves  <- as.data.frame(cbind(resumenP, resumenQ, resumenR, resumenS, resumenT))
resumenWaves

# write.xlsx(resumenWaves, "AnalisisInfFMM3D-M3-500cobWaves.xlsx", colNames = TRUE, rowNames = TRUE,
#            colWidths = 10)



