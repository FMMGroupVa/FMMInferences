
library(R.utils)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rgl)
source("numericInferences.R")
source("MatrizPrecision.R")
source("simulationFunctions.R")
sourceDirectory("funcionesFMM/")

################################################################################

nExperiments <- 1000
nCores <- 64
M <- 0; A <- 1; alpha <- 5;
sigmaValues <- c(0.025, 0.05)
nValues <- c(51, 101, 201, 401)
#nValues <- c(50, 100)
betaValues <- c(pi, 3*pi/2)
omegaValues <- c(0.05, 0.1, 0.15, 0.2)
#omegaValues <- c(0.05, 0.1, 0.15)

parameterGrid <- expand.grid(sigmaValues, nValues, betaValues, omegaValues)
#parameterGrid2 <- expand.grid(sigmaValues, nValues, betaValues, omegaValues2)

groupOfExperiments <- function(arg){
  arg <- as.numeric(arg)
  sigmaVal <- arg[1]; nVal <- arg[2]; betaVal <- arg[3]; omegaVal <- arg[4]
  experiment <- covertureExperimentFMM1(M = 0, A = 1, alpha = 5, beta = betaVal,
                                        omega = omegaVal, sigma = sigmaVal, n = nVal,
                                        nExperiments = nExperiments, conf = 0.95)
  return(experiment)
}

X <- asplit(parameterGrid, 1)

red <- parallel::mclapply(X = X, FUN = groupOfExperiments, mc.cores = nCores)

red2 <- bind_rows(red, .id = "column_label")


write.table(red2, "FMM1-AllExperiments.csv", row.names = F)

################################################################################


simulFMM1 <- red2

coverages <- simulFMM1 %>% group_by(sigma, n, beta, omega) %>% 
  summarise(alphaF0 = sum(alphaIn, na.rm = T)*100/nExperiments, 
            omegaF0 = sum(omegaIn, na.rm = T)*100/nExperiments, 
            alphaHess = sum(alphaInHess, na.rm = T)*100/nExperiments, 
            omegaHess = sum(omegaInHess, na.rm = T)*100/nExperiments,
            alphaHP = sum(alphaInProf, na.rm = T)*100/nExperiments, 
            omegaHP = sum(omegaInProf, na.rm = T)*100/nExperiments,
            lenAlphaF0 = 4*mean(sdalphaF0),
            lenOmegaF0 = 4*mean(sdomegaF0),
            lenAlphaHess = 4*mean(sdHessalpha),
            lenOmegaHess = 4*mean(sdHessomega),
            lenAlphaProf = 4*mean(sdProfalpha),
            lenOmegaProf = 4*mean(sdProfomega))
coverages
# write.xlsx(coverages, "coveragesFMM1.xlsx", colNames = TRUE, rowNames = FALSE, 
#            colWidths = 8)






covsAlphas <- coverages[,c(1:4, 5, 7, 9)]
covsAlphas <- covsAlphas %>% pivot_longer(alphaF0:alphaHP, names_to = "Estimator", 
                                          values_to = "Mean_Coverage") 
covsAlphas$beta <- as.character(round(covsAlphas$beta, 2))
covsAlphas$beta <- factor(covsAlphas$beta, labels = c("pi", "pi/2"))
covsAlphas$Estimator <- factor(covsAlphas$Estimator)

covsAlphas %>%
  ggplot(aes(x = n, y = Mean_Coverage, groups = Estimator, colour = Estimator)) +
  geom_hline(yintercept = 95, size = 1.4, alpha = 0.7, col = "grey") +
  geom_line() +
  geom_point() +
  scale_x_continuous(trans='log2', breaks=c(50, 100, 200, 400)) +
  
  xlab("Sample size (n)") + ylab("Mean coverage (%)") +
  facet_grid(sigma + beta ~ omega, margins = "VS") +
  theme_light() +
  guides(colour = guide_legend(nrow = 1)) +
  theme(axis.title = element_text(face = "bold"), legend.position = "bottom", 
        legend.title = element_blank(), legend.text = element_text(face = "bold"))
