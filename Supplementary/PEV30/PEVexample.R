
# Patient 5: PEV30 example. ICs for peak latency

# IMPORTANT: set the working directory to this file location 
# setwd(".../FMMInferences-main")

library(readxl)
source("auxMultiFMM.R")

dataRE <- read.csv("RightEyePEVda30.csv")
dataLE <- read.csv("LeftEyePEVda30.csv")
signalDuration <- 375

paramsPerSignalRE<-fitMultiFMM(vDataMatrix=dataRE, nBack=4, maxIter=100, weightError = TRUE, showPredeterminedPlot = T)
paramsPerSignalLE<-fitMultiFMM(vDataMatrix=dataLE, nBack=4, maxIter=100, weightError = TRUE, showPredeterminedPlot = T)

CI_RE <- confint(paramsPerSignal = paramsPerSignalRE, mData = dataRE)
CI_LE <- confint(paramsPerSignal = paramsPerSignalLE, mData = dataLE)
CIlist <-list(CI_RE, CI_LE)

n <- 5000 # Bootstrap samples
peak1<-array(0,c(4,n))
peak2<-array(0,c(4,n))
sdi = function(ic) unname((ic[3]-ic[2])/2)

# Generation of bootstrap samples: 
for (s in 1:2) {
  
  CIS <- CIlist[[s]]
  alphaBoot <- rnorm(n = n, mean = CIS[2,"alpha_1"], sd = sdi(CIS[,"alpha_1"]))
  omegaBoot <- rnorm(n = n, mean = CIS[2,"omega_1"], sd = sdi(CIS[,"omega_1"]))
  delta1Boot <- rnorm(n = n, mean = CIS[2,"delta_1_Ch1"], sd = sdi(CIS[,"delta_1_Ch1"]))
  gamma1Boot <- rnorm(n = n, mean = CIS[2,"gamma_1_Ch1"], sd = sdi(CIS[,"gamma_1_Ch1"]))
  delta2Boot<- rnorm(n = n, mean = CIS[2,"delta_1_Ch2"], sd = sdi(CIS[,"delta_1_Ch2"]))
  gamma2Boot<- rnorm(n = n, mean = CIS[2,"gamma_1_Ch2"], sd = sdi(CIS[,"gamma_1_Ch2"]))

  betas_Ch1<-(atan2(-gamma1Boot,delta1Boot))%%(2*pi)
  betas_Ch2<-(atan2(-gamma2Boot,delta2Boot))%%(2*pi)
  
  # Peak locations (prominent wave)
  peak1[s,]<-(alphaBoot+2*atan2((1/omegaBoot)*sin(-betas_Ch1/2), cos(-betas_Ch1/2)))%%(2*pi)*signalDuration/(2*pi)
  peak2[s,]<-(alphaBoot+2*atan2((1/omegaBoot)*sin(-betas_Ch2/2), cos(-betas_Ch2/2)))%%(2*pi)*signalDuration/(2*pi)
}

# Bootstrap confidence intervals for peak latency
# Right eye
REPeak1CI <- round(quantile(peak1[1,],probs=c(0.025,0.5,0.975)),1) # 113[107.5, 117.9]
REPeak2CI <- round(quantile(peak2[1,],probs=c(0.025,0.5,0.975)),1) # 112[106.7, 116.3]
# Left eye
LEPeak1CI <- round(quantile(peak1[2,],probs=c(0.025,0.5,0.975)),1) # 122[120.7,124.2] 
LEPeak2CI <- round(quantile(peak2[2,],probs=c(0.025,0.5,0.975)),1) # 121[119.1, 122.4]
# Difference
DiffPeak1CI <- round(quantile((peak1[2,]-peak1[1,]),probs=c(0.025,0.5,0.975)),1)# [4.3, 15.3]  
DiffPeak2CI <- round(quantile((peak2[2,]-peak2[1,]),probs=c(0.025,0.5,0.975)),1) # [4.1, 14.3]  

cat("Right eye peaks CIs:\n\t 5%\t 50%\t 95%",
    "\nCh1:\t", REPeak1CI[1],"\t", REPeak1CI[2],"\t", REPeak1CI[3],"\t",
    "\nCh2:\t", REPeak2CI[1],"\t", REPeak2CI[2],"\t", REPeak2CI[3],"\t") 

cat("Left eye peaks CIs:\n\t 5%\t 50%\t 95%",
    "\nCh1:\t", LEPeak1CI[1],"\t", LEPeak1CI[2],"\t", REPeak1CI[3],"\t",
    "\nCh2:\t", LEPeak2CI[1],"\t", LEPeak2CI[2],"\t", REPeak2CI[3],"\t") 

# Right eye
cat("Right-Left peak difference CIs:\n\t 5%\t 50%\t 95%",
    "\nCh1:\t", DiffPeak1CI[1],"\t", DiffPeak1CI[2],"\t", DiffPeak1CI[3],"\t",
    "\nCh2:\t", DiffPeak2CI[1],"\t", DiffPeak2CI[2],"\t", DiffPeak2CI[3],"\t") 



  

