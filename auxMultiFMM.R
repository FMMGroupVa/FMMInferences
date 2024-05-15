#### Dependencies ####
require("RColorBrewer")
require(R.utils)
sourceDirectory("funcionesFMM/")
source("auxMultiFMMPlot.R")
source("precissionMatrix.R")

#### Fit of multiFMM model functions ####
fitMultiFMM <- function(vDataMatrix, nBack = 5, maxIter = 10, weightError = TRUE,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMin = 0.001, omegaMax = 1,
                      omegaGrid = exp(seq(log(max(omegaMin, omegaMin)), log(1), length.out = lengthOmegaGrid)),
                      parallelize = TRUE, confidenceLevel = 0.95, plotToFile = F, filename = NA){

  nSignals <- ncol(vDataMatrix)
  nObs <- nrow(vDataMatrix)
  #### Preparations for the fit ####
  gridList <- precalculateBase(alphaGrid = alphaGrid, omegaGrid = omegaGrid,
                               timePoints = seqTimes(nObs))

  # The results are stored in:
  paramsPerWave <- replicate(nSignals, simplify = FALSE,
                             data.frame(M=double(), A=double(), Alpha=double(), Beta=double(), Omega=double()))
  fittedWaves<-lapply(1:nSignals,function(x){
    vData<-vDataMatrix[,x]; vData<-vData[!is.na(vData)];
    return(array(0,c(length(vData),nBack)))})

  # The MSE must be balanced in the fit of multiple signals
  errorWeights<-rep(1,nSignals); totalMSE<-rep(Inf, maxIter)

  # Backfitting algorithm: iteration
  currentBack<-1; continueBackfitting<-TRUE
  cat(paste("Backfitting: ",sep=""))

  while(continueBackfitting) {
    cat(paste(currentBack," ",sep=""))
    fmmObjectList<-list()
    paramsPerWave <- replicate(nSignals, simplify = FALSE,
                               data.frame(M=double(), A=double(),
                                          Alpha=double(), Beta=double(), Omega=double()))
    # Backfitting algorithm: component
    for(j in 1:nBack){
      #### First Step: determine optimal common parameters (alpha and omega) ####
      optimalParams <- optimizeAlphaOmega(vDataMatrix = vDataMatrix, baseGrid = gridList,
                                          fittedWaves = fittedWaves, currentComp = j,
                                          omegaMax = omegaMax, errorWeights = errorWeights)

      #### Second Step: fit single FMM wave in each signal with common parameters ####
      for(signalIndex in 1:nSignals){
        vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
        vData<-vData - apply(as.matrix(fittedWaves[[signalIndex]][,-j]), 1, sum)
        fmmObjectList[[signalIndex]]<-optimizeOtherParameters(vData = vData, fixedParams = optimalParams)

        fittedWaves[[signalIndex]][,j]<-getFittedValues(fmmObjectList[[signalIndex]])
      }

      #### Get FMM parameters per wave ####
      params<-data.frame("M"=sapply(fmmObjectList,getM),"A"=sapply(fmmObjectList,getA),"Alpha"=sapply(fmmObjectList,getAlpha),
                         "Beta"=sapply(fmmObjectList,getBeta),"Omega"=sapply(fmmObjectList,getOmega), "Var"=rep(NA,nSignals))
      paramsPerWave<-lapply(1:nSignals, function(x) rbind(paramsPerWave[[x]], params[x,]))
      paramsPerWave<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], paramsPerSignal=paramsPerWave[[x]]))

      #### Error is weighted across signals ####
      sigma<-sapply(1:nSignals, function(x){sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]]), 1, sum))^2)/(nObs-1)})
      if(weightError) errorWeights<-1/sigma
      else errorWeights<-rep(1,nSignals)
    }
    totalMSE[currentBack]<-sum(sapply(1:nSignals, function(x){sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]]), 1, sum))^2)}))

    # Check if backfitting should stop
    stopCondition1<-currentBack>=maxIter
    stopCondition2<-ifelse(currentBack>1, totalMSE[currentBack]>totalMSE[currentBack-1], FALSE)
    if(stopCondition1 | stopCondition2){
      stopCriteria<-ifelse(stopCondition1, "Maximum Iterations", "Minimum MSE")
      cat(paste("  Stop: ", stopCriteria, "\n", sep=""))
      continueBackfitting<-FALSE
    }else{
      currentBack<-currentBack+1
    }
  }

  # Reorder waves by mean Var across channels
  meanVar <- apply(matrix(unlist(lapply(paramsPerWave, FUN = function(x){
    return(x$Var)
  })), ncol = nSignals, byrow = F), 1, mean)
  
  waveOrder <- order(meanVar, decreasing = T)

  for (i in 1:nSignals) {
    paramsPerWave[[i]] <- (paramsPerWave[[i]])[waveOrder,]
    fittedWaves[[i]] <- (fittedWaves[[i]])[,waveOrder]
  }
  
  
  plotMultiFMM(vDatai = vDataMatrix, fittedWaves = fittedWaves, currentBack = currentBack,
               leadNames = colnames(vDataMatrix), paramsPerSignal = paramsPerWave,
               plotToFile = plotToFile, filename = filename)

  # Unname waves and stop parallelized cluster
  for(i in 1:nSignals) rownames(paramsPerWave[[i]])<-1:nBack

  # Confidence Intervals calculus
  CIs <- confint(paramsPerSignal = paramsPerWave, mData = vDataMatrix,
                 nBack = nBack, nSignals = nSignals,
                 compNames = 1:nBack,  confidenceLevel = 0.95)
  #### Return results ####
  return(list(paramsPerWave = paramsPerWave, Confints = CIs))
}
#### Internal multiFMM functions ####
## MultiFMM, first step: optimize common parameters

# logLik3DFMM1 <- function(alphaOmegaParameters, vData, timePoints, sigmas = rep(1, ncol(vData))) {
#   alpha <- as.numeric(alphaOmegaParameters[1])
#   omega <- as.numeric(alphaOmegaParameters[2])
#   tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
#   fittedValues <- apply(vData, 2, function(x){
#     dM <- cbind(rep(1, length(x)), cos(tStar), sin(tStar))
#     mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%x
#     yFit <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
#     return(yFit)
#   })
#   residuales <- (vData - fittedValues)^2
#   logL <- -sum(t(sigmas^2*t(residuales))/2) # Cada columna por su sigma^2
#   return(logL)
# }

step2_commonAlphaOmega <- function(initialParams, vDataMatrix, weights){

  initialAlpha<-as.numeric(initialParams[1])
  initialOmega<-as.numeric(initialParams[2])
  timePoints <- seqTimes(length(vDataMatrix[,1]))

  tStar <- 2*atan(initialOmega*tan((timePoints-initialAlpha)/2))

  # Calculates estimated component per lead
  fittedValues <- apply(vDataMatrix, 2, function(x, tStar){
    dM <- cbind(rep(1, length(x)), cos(tStar), sin(tStar))
    mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%x
    yFit <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
    return(yFit)
  }, tStar = tStar)

  # Weighted residuals
  residuales <- (vDataMatrix - fittedValues)^2
  RSS <- sum(t(weights*t(residuales)))

  return(RSS)
}

optimizeAlphaOmega<-function(vDataMatrix, baseGrid, fittedWaves, currentComp,
                             errorWeights, omegaMax = 0.7){

  residualsMatrix<-as.matrix(vDataMatrix)
  for(signalIndex in 1:ncol(vDataMatrix)){
    vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
    residualsMatrix[,signalIndex] <- vData - apply(as.matrix(fittedWaves[[signalIndex]][,-currentComp]), 1, sum)
  }
  # Grid step. RSS is a weighted mean of the RSS(i), i in cols(vDataMatrix)
  step1 <- lapply(FUN = step1FMM3D, X = baseGrid, vDataMatrix = residualsMatrix,
                  weights = errorWeights)

  step1 <- matrix(unlist(step1), ncol=3, byrow=T)

  bestParamsIndex <- which.min(step1[,3])

  alpha <- step1[bestParamsIndex, 1]
  omega <- step1[bestParamsIndex, 2]

  # Post-optimization. Depends on (alpha, omega)
  nelderMead <- optim(par = c(alpha, omega), fn = step2_commonAlphaOmega,
                      vDataMatrix = residualsMatrix, method = "L-BFGS-B",
                      lower = c(-2*pi, 0.0001), upper = c(4*pi, omegaMax), weights = errorWeights)

  return(nelderMead$par[1:2])
}

## MultiFMM, second step: determine non-common parameters
optimizeOtherParameters <- function(vData, fixedParams, timePoints = seqTimes(length(vData))){

  alpha <- as.numeric(fixedParams[1])
  omega <- as.numeric(fixedParams[2])
  tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))

  dM <- cbind(rep(1, length(tStar)), cos(tStar), sin(tStar))
  mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%vData
  M <- mDeltaGamma[1]
  delta <- mDeltaGamma[2]
  gamma <- mDeltaGamma[3]

  fittedFMMvalues <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
  SSE <- sum((vData-fittedFMMvalues)^2)

  fmmObject<-FMM(
    M = mDeltaGamma[1],
    A = sqrt(delta^2+gamma^2),
    alpha = alpha,
    beta = atan2(-gamma, delta) + alpha,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = PV(vData, fittedFMMvalues),
    nIter = 0
  )

  fmmObject@nPeriods<-1
  fmmObject@data<-vData

  return(fmmObject)
}

step2FMM_fixedAlphaOmega <- function(parameters, vData, timePoints, fixedAlpha, fixedOmega){
  nObs <- length(timePoints)
  modelFMM <- parameters[1] + parameters[2] *
    cos(parameters[4]+2*atan2(fixedOmega*sin((timePoints - fixedAlpha)/2),
                              cos((timePoints - fixedAlpha)/2)))
  residualSS <- sum((modelFMM - vData)^2)/nObs
  rest3 <- parameters[2] > 0  # A > 0
  if(rest3)
    return(residualSS)
  else
    return(Inf)
}

#### Other multiFMM useful functions ####
recalculateMA<-function(vDatai, paramsPerSignal){

  nObs<-length(vDatai)

  ## Extract Params. and calculate cosPhi
  alpha<-paramsPerSignal$Alpha; beta<-paramsPerSignal$Beta
  omega<-paramsPerSignal$Omega; maxComp<-length(alpha)
  cosPhi<-t(calculateCosPhi(alpha = alpha, beta = beta,
                            omega = omega, timePoints = seqTimes(nObs)))

  ## Recalculate M and As
  currentFormula<-as.formula(paste("vDatai~", paste("cosPhi[",1:maxComp,",]",sep="",collapse="+")))
  regression<-lm(currentFormula)
  recalculatedM<-regression$coefficients[1]
  recalculatedAs<-as.vector(regression$coefficients[-1])

  ## Calculate R2, including other waves
  individualR2<-PVj(vData = vDatai, timePoints = seqTimes(nObs),
                    alpha = paramsPerSignal$Alpha,
                    beta = paramsPerSignal$Beta, omega = paramsPerSignal$Omega)
  paramsPerSignal$Var<-individualR2

  mResults<-rep(NA,maxComp); aResults<-rep(NA,maxComp)
  mResults<-rep(recalculatedM, maxComp); aResults<-recalculatedAs
  paramsPerSignal$M<-mResults; paramsPerSignal$A<-aResults

  ## Negative A => Translation on the beta parameter
  aValidContion <- aResults>0 | is.na(aResults)
  if(any(!aValidContion)){
    paramsPerSignal$Beta[!aValidContion]<-(paramsPerSignal$Beta[!aValidContion]+pi)%%(2*pi)
    paramsPerSignal$A[!aValidContion]<- -paramsPerSignal$A[!aValidContion]
  }

  return(paramsPerSignal)
}

confint <- function(paramsPerSignal, mData,  nBack, nSignals,
                    confidenceLevel = 0.95, compNames = 1:nBack){

  estimatesArray <- array(0, c(nBack, 5, nSignals))
  for (k in 1:nSignals){ estimatesArray[,,k] <- as.matrix(paramsPerSignal[[k]][,-6]) }
  alphasHat <- estimatesArray[,3,1]; omegasHat <- estimatesArray[,5,1]
  AsHat <- estimatesArray[,2,]; betasHat <- estimatesArray[,4,]
  deltasHat <- c(AsHat*cos(betasHat)); gammasHat <- c(-AsHat*sin(betasHat))
  estimatedParameters <- c(deltasHat, gammasHat, alphasHat, omegasHat)

  # Sigma estimation
  r2 <- sapply(1:length(paramsPerSignal), function(x) sum(paramsPerSignal[[x]][substr(rownames(paramsPerSignal[[x]]),1,1)!="X",]$Var))
  sigmaHat <- sqrt(sum((1-r2)*apply(mData, 2, function(vData){sum((vData - mean(vData))^2)}))/(nSignals*nrow(mData) - (1+(1+nSignals)*2*nBack)))

  # CIs with F0
  tF0F0 <- precissionMatrixFMM3d(estimatesArray, seqTimes(nrow(mData)))
  seHat <- sqrt(diag(sigmaHat^2*solve(tF0F0)))[-1]
  lInf <- estimatedParameters - qnorm(0.5+confidenceLevel/2)*seHat
  lSup <- estimatedParameters + qnorm(0.5+confidenceLevel/2)*seHat

  # Format (row and column names)
  leadNames <- colnames(mData)
  namesDeltas <- paste0("delta_",rep(compNames, nSignals), "_",rep(leadNames, each = nBack))
  namesGammas <- paste0("gamma_",rep(compNames, nSignals), "_",rep(leadNames, each = nBack))
  namesAlphas <- paste0("alpha_",compNames)
  namesomega <- paste0("omega_",compNames)
  fullNames <- c(namesDeltas, namesGammas, namesAlphas, namesomega)

  CIs <- rbind(lInf, estimatedParameters, lSup)
  colnames(CIs) <- fullNames
  rownames(CIs) <- c(paste0("Lower (", 100*confidenceLevel, "%)"), "Estimate", paste0("Upper (", 100*confidenceLevel, "%)"))
  return(CIs)
}
