
# IMPORTANT: set the working directory to this file location 
# setwd("---/FMMInferences-main")

# Libraries for source code reading
source("auxMultiFMM.R") # Codes for 3DFMM estimation + CIs

# Simulated data with common sigma
exampleData <- read.csv("exampleData.csv")

################################################################################

# 1. Estimation of the parameters with -fitMultiFMM- function.

# Arguments: 
#     vDataMatrix - Data matrix of dimenstion num observations x num Channels
#     nBack - Number of waves 
#     maxIter - Number of iterations in the backfitting algorithm. To ensure the
#               convergence, a high number of iterations are needed.
#     showPredeterminedPlot - Plot with prediction + waves 

paramsPerSignal <- fitMultiFMM(vDataMatrix = exampleData, nBack = 5, maxIter = 10,
                               plotToFile = F, showPredeterminedPlot = T)

################################################################################

# 2. Confidence intervales for FMM parameters

# Arguments: 
#     paramsPerSignal - output of fitMultiFMM function
#     mData - data matrix fitted. Used for the estimate of sigma
#     confidenceLevel - confidence level for the CIs

CIs <- confintFMM(paramsPerSignal = paramsPerSignal, mData = exampleData,
                  confidenceLevel = 0.95)

print(paramsPerSignal[[2]]) # Estimated parameters (Lead II)
print(round(CIs[,81:90],3)) # CIs for alphas and omegas

################################################################################

# 3. To show a selection of channels use the function plotMultiFMM(...)

FMM3D_output <- fitMultiFMM(vDataMatrix = exampleData, nBack = 5, maxIter = 10,
                            plotToFile = F, showPredeterminedPlot = F) # Do not show predetermined plot

# User must give data matrix, parameters output, and specify what channels 
# are required and how to plot them:

# Fitted signal plots
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output, 
             plotToFile = F, nPlotCols = 3, channels = 1:6)

# Component plots
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output, 
             plotToFile = F, nPlotCols = 3, channels = 1:6, components = T)


