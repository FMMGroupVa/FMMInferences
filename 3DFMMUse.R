
# Libraries for source code reading
require(R.utils)
source("auxMultiFMM.R") # Codes for 3DFMM estimation + CIs
sourceDirectory("funcionesFMM/")

# Simulated data with common sigma
exampleData <- read.csv("exampleData.csv")

FMM3D_output <- fitMultiFMM(vDataMatrix = exampleData, nBack = 5, maxIter = 5,
                            parallelize = TRUE, confidenceLevel = 0.95, plotToFile = F,
                            showPredeterminedPlot = T)

paramsPerSignal <- FMM3D_output[[1]]; CIs <- FMM3D_output[[2]] # Estimates / CIs

print(paramsPerSignal[[2]]) # Estimated parameters (Lead II)
print(round(CIs[,81:90],3)) # CIs for alphas and omegas

################################################################################

# To show a selection of channels use the functions

FMM3D_output <- fitMultiFMM(vDataMatrix = exampleData, nBack = 5, maxIter = 5,
                            parallelize = TRUE, confidenceLevel = 0.95, 
                            plotToFile = F, showPredeterminedPlot = F) # Do not show predetermined plot

# User must give data matrix, parameters output, and specify what channels 
# are required and how to plot them:

# Fitted signal plots
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output[[1]], 
             plotToFile = F, nPlotCols = 2, channels = 1:4)
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output[[1]], 
             plotToFile = T, nPlotCols = 2, channels = 1:4, 
             filename = "ex1.png") # Plotted to file ex1.png

# Component plots
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output[[1]], 
             plotToFile = F, nPlotCols = 2, channels = 1:4, components = T)
plotMultiFMM(vData = exampleData, paramsPerSignal = FMM3D_output[[1]], 
             plotToFile = T, nPlotCols = 2, channels = 1:4, components = T, 
             filename = "ex2.png") # Plotted to file ex2.png


