
# Libraries for source code reading
require(R.utils)
source("auxMultiFMM.R") # Codes for 3DFMM estimation + CIs
sourceDirectory("funcionesFMM/")

# Simulated data with common sigma
exampleData <- read.csv("exampleData.csv")

FMM3D_output <- fitMultiFMM(vDataMatrix = exampleData, nBack = 5, maxIter = 15, 
                            parallelize = TRUE, confidenceLevel = 0.95)

paramsPerSignal <- FMM3D_output[[1]]; CIs <- FMM3D_output[[2]] # Estimates / CIs

print(paramsPerSignal[[2]]) # Estimated parameters (Lead II)
print(round(CIs[,81:90],3)) # CIs for alphas and omegas

