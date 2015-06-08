TableCompare <- function(mvals, nreps, nimFilt, otherFilt, otherPack, filtType, like=T){
  timelist <- data.frame()
  nimTempTimes <- rep(NA, nreps)
  nimTempLikes <- rep(NA, nreps)
  othTempTimes <- rep(NA, nreps)
  othTempLikes <- rep(NA, nreps)
  for(m in mvals){
    for(n in 1:nreps){
      nimTempTimes[n] <- sum(system.time(nimTempLikes[n] <- nimFilt$run(m))[1:2])
      if(like == T){
        othTempTimes[n] <- sum(system.time(othTempLikes[n] <- otherFilt(m))[1:2])
      }
      else{
        othTempTimes[n] <- sum(system.time(otherFilt(m))[1:2])
      }
    }    
    
    if(like == T){
      timelist <- rbind(timelist,cbind("Method" = c(paste("NIMBLE",filtType, sep=" "),
                                                    paste(otherPack,filtType,sep=" ")),
                                       "N. Particles" =c(m,m),
                                       "Mean Log Lik." = round(c(mean(nimTempLikes),  mean(othTempLikes)),2),
                                       rbind(summary(nimTempTimes), summary(othTempTimes))[,c(1,3,4,6)]))
    }
    else{
      timelist <- rbind(timelist,cbind("Method" = c(paste("NIMBLE",filtType, sep=" "),
                                                    paste(otherPack,filtType,sep=" ")),
                                       "N. Particles" =c(m,m),
                                       rbind(summary(nimTempTimes), summary(othTempTimes))[,c(1,3,4,6)]))
    }
  }
  return(timelist)
}

PlotCompare <- function(plotm, nimFilt, otherOut, otherPack, filtType, othWts, kalmanDat=NA){
  plotMat <- data.frame(value = c(as.matrix(nimFilt$mv, 'x'),
                                  otherOut),
                        package=c(rep("Nimble",plotm),rep(otherPack,plotm)),
                        wts= c(as.matrix(nimFilt$mv, 'wts'), othWts))
  
  #think about what binwidth should be
  compPlot <- ggplot(plotMat, aes(x=value))+
    geom_histogram(aes(y=..density.., weight=wts),
                   colour = "darkgreen", fill = "white", binwidth = 0.02)+ 
    facet_wrap(~package, nrow=2, as.table=F)+ 
    ggtitle(paste(filtType, "Posteriors for x",sep=" "))
  
  if(!any(is.na(kalmanDat)))
    compPlot <- compPlot + geom_line(data=kalmanDat, aes(x=nums, y=dens))

  return(suppressWarnings(compPlot))
}
