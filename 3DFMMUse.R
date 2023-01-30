
library(openxlsx)
require(R.utils)
source("auxMultiFMM.R")
sourceDirectory("funcionesFMM/")

normBeat <- read.xlsx("param_config.xlsx", "NORM", )
leadNames <- normBeat[,1]
leadParams <- normBeat[,-1]
m <- 5
n <- 200
d <- length(leadNames)
commonSigma <- 45

paramArray <- array(0, c(m, 5, d))
for (k in 1:d){ paramArray[,,k] <- matrix(as.numeric(leadParams[k,1:(5*m)]), ncol = 5, byrow = T) }

simulatedData <- apply(paramArray, 3, function(x){generateFMM(x[1,1], x[,2], x[,3], x[,4], x[,5], 
                                                              sigmaNoise = commonSigma, length.out = n, plot = F)$y})
colnames(simulatedData) <- leadNames
FMM3D_output <- fitMultiFMM(vDataMatrix = simulatedData, nBack = 5, maxIter = 10)

paramsPerSignal <- FMM3D_output[[1]]; CIs <- FMM3D_output[[2]]
paramsPerSignal
print(round(CIs[,81:90],3))




