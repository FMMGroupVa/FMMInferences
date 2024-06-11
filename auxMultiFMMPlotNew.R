
#### Plot multiFMM functions ####
emptyPlotFun<-function() plot(1, type = "n", axes=FALSE, xlab="", ylab="")

plotMultiFMM<-function(vDatai, timePoints = seqTimes(nrow(vDatai)), 
                       paramsPerSignal, channels = 1:ncol(vDatai), components = F, nPlotCols = 5, 
                       filename=NA, leadNames=1:length(channels), 
                       path="./", plotToFile=TRUE){

  vDatai <- vDatai[,channels]
  paramsPerSignal <- paramsPerSignal[channels]
  
  nObs <- nrow(vDatai)
  nSignals <- ncol(vDatai)
  leadNames <- colnames(vDatai)
  maxComp <- nrow(paramsPerSignal[[1]])
  nPlotRows = ceiling(nSignals/nPlotCols) # Number of row for fixed columns
  if(plotToFile){
    if(!is.na(filename)){
      png(filename=paste(path,filename,".png",sep=""), type = "cairo-png",
          width=ifelse(nSignals<=12 | nSignals%%3!=0,
                       200*nSignals,50*nSignals),
          height=ifelse(nSignals<=12, 600,
                        ifelse(nSignals%%3!=0, 600*round(nSignals/12), 700*round(nSignals/12))))
    }else{
      stop("Filename must be provided if plotToFile=TRUE")
    }
  }

  fittedWaves <- list() # Calcularlas
  
  for(s in 1:nSignals){

    fittedWaves[[s]] <- matrix(0, nrow = nObs, ncol = maxComp)
    for (k in 1:maxComp) {
      pars <- as.numeric((paramsPerSignal[[s]])[k,])
      fittedWaves[[s]][,k] <- pars[2]*cos(pars[4] + 2*atan(pars[5]*tan((timePoints-pars[3])/2)))
    }
  }

  #predictedSignals <-  lapply(fittedWaves, function(x){apply(x, 2, sum)})

  if(components){
    layout(matrix(c(1:(nPlotRows * nPlotCols), rep(nPlotRows*nPlotCols+1,nPlotCols)), 
                  nrow=nPlotRows+1, byrow = T),
           heights=c(rep(3, times = nPlotRows),1))
    
    for (k in 1:(nPlotRows* nPlotCols)) {
      if(k<=nSignals){
        plotMultiFMM_Comps(vDatai=vDatai[,k], fittedWaves = fittedWaves[[k]],
                           paramsPerSignal = paramsPerSignal[[k]],
                           filename=filename, leadName=leadNames[k],
                           plotLegend=FALSE)
      }else{
        emptyPlotFun()
      }
    }
    emptyPlotFun()
    
    usedColors<-c(brewer.pal(n = 9, name = "Set1"), "aquamarine", "saddlebrown")
    usedColors<-usedColors[1:nrow(paramsPerSignal[[1]])]
    names(usedColors)<-1:nrow(paramsPerSignal[[1]])
    legend(x = "top",inset = 0, legend = names(usedColors), col=usedColors,
           lwd=5, cex=1.2, horiz = TRUE)

  }else{
    par(mfrow = c(nPlotRows, nPlotCols))
    for (k in 1:nSignals) {
      plotMultiFMM_Sum(vDatai=vDatai[,k], fittedWaves = fittedWaves[[k]],
                       paramsPerSignal = paramsPerSignal[[k]],
                       filename=filename, leadName=leadNames[k])
    }
  }
  par(mfrow = c(1, 1))

  if(plotToFile) dev.off()
}

plotMultiFMM_Sum<-function(vDatai, fittedWaves, paramsPerSignal, leadName,
                           filename=NA, path="./", plotToFile=FALSE){
  par(mar = c(2,2,3,1))
  
  nObs<-length(vDatai)
  assignedCondition<-rep(TRUE,nrow(paramsPerSignal))
  
  assignedResults<-paramsPerSignal[assignedCondition,]
  assignedWavesSum<-generateFMM(M = assignedResults$M[1], A = assignedResults$A,
                                alpha = assignedResults$Alpha, beta = assignedResults$Beta,
                                omega = assignedResults$Omega, length.out = nObs, plot = F)$y
  
  ## Colors definition: color linked to assigned wave
  totalR2<-sum(assignedResults$Var)
  
  if(plotToFile){
    png(filename=paste(path,"/02 Results/Plots/",filename,".png",sep=""),
        width=1200, height=900)
  }
  
  mainText<-leadName
  
  assignedWavesSumPlot<-assignedWavesSum
  
  yLimits<-c(min(assignedWavesSumPlot,vDatai), max(assignedWavesSumPlot,vDatai))
  plot(1:nObs,vDatai,type="l",ylim=yLimits,cex.main=1.6,
       main=mainText)
  lines(1:nObs,assignedWavesSumPlot,col=4,lwd=2)
  legend("bottomright",legend=sprintf(totalR2*100, fmt = "R2=%#.1f %%"), col=NA, lty=1,
         lwd=NA, pch=NA, text.width=40, x.intersp=0.25, bty = "n", cex=1.5,
         inset=c(0.2,0))
  
  if(plotToFile) dev.off()
}

plotMultiFMM_Comps<-function(vDatai, fittedWaves, paramsPerSignal,
                             leadName, filename=NA, path="./",
                             plotLegend=TRUE, plotToFile=FALSE){
  par(mar = c(2,2,1,1))
  
  maxComp<-nrow(paramsPerSignal)
  assignedCondition<-rep(TRUE,nrow(paramsPerSignal))
  
  ## Colors definition: color linked to assigned wave
  assignedWavesColors<-c(brewer.pal(n = 9, name = "Set1"), "aquamarine", "saddlebrown")
  assignedWavesColors<-assignedWavesColors[1:nrow(paramsPerSignal)]
  otherWavesColors<-rev(brewer.pal(n = 7, name = "Set2"))
  
  if(plotToFile){
    png(filename=paste(path,"/02 Results/Plots/",filename,".png",sep=""),
        width=1200, height=900)
  }
  
  ## Components plot
  yLimits<-c(min(sapply(1:maxComp, function(x) fittedWaves[,x]-median(fittedWaves[,x]))),
             max(sapply(1:maxComp, function(x) fittedWaves[,x]-median(fittedWaves[,x]))))
  currentColor<-ifelse(assignedCondition[1], assignedWavesColors[1], otherWavesColors[1])
  
  plot(fittedWaves[,1]-median(fittedWaves[,1]), lty=ifelse(assignedCondition[1],1,2),
       col=currentColor, type="l", lwd=3, ylim=yLimits)
  if(assignedCondition[1]) usedColors<-c(currentColor)
  
  if(maxComp>1){
    notAssignedColorIndex<-ifelse(!assignedCondition[1],2,1)
    assignedColorIndex<-ifelse(assignedCondition[1],2,1)
    for(i in 2:maxComp){
      currentColor<-ifelse(assignedCondition[i], assignedWavesColors[assignedColorIndex],otherWavesColors[notAssignedColorIndex])
      if(assignedCondition[i]){
        if(!exists("usedColors")){usedColors<-c(currentColor)
        }else{usedColors<-c(usedColors, currentColor)}
        assignedColorIndex<-assignedColorIndex+1
      }else{
        notAssignedColorIndex<-notAssignedColorIndex+1
      }
      lines(fittedWaves[,i]-median(fittedWaves[,i]),col=currentColor,
            lwd=3, lty=ifelse(assignedCondition[i],1,2))
    }
  }
  
  ## Legend of the components plot
  if(plotLegend){
    legend("bottomright",legend=paste("Wave",1:maxComp),col=usedColors,lty=rep(1,maxComp),
           lwd=rep(2,maxComp),pch=rep(NA,maxComp), cex = 1, text.width=25, x.intersp=0.25,
           y.intersp = 0.1, bty = "n")
  }
  
  if(plotToFile) dev.off()
}




