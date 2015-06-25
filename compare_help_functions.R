TableCompare <- function(mvals, nreps, nimFilt, otherFilt, otherPack, filtType, like=T){
  timelist <- data.frame()
  nimTempTimes <- rep(NA, nreps)
  nimTempLikes <- rep(NA, nreps)
  othTempTimes <- rep(NA, nreps)
  othTempLikes <- rep(NA, nreps)
  for(m in mvals){
    for(n in 1:nreps){
      if(like == T){
        nimTempTimes[n] <- sum(system.time(nimTempLikes[n] <- nimFilt$run(m))[1:2])
     
        othTempTimes[n] <- sum(system.time(othTempLikes[n] <- otherFilt(m))[1:2])
      }
      else{
        othTempTimes[n] <- sum(system.time(otherFilt(m))[1:2])
        nimTempTimes[n] <- sum(system.time(nimFilt$run(m))[1:2])
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


PMCMCTableCompare <- function(rvals,mvals, nreps, otherFilt, otherPack, filtType){
  timelist <- data.frame()
  nimTempTimes <- rep(NA, nreps)
  nimTempLikes <- rep(NA, nreps)
  othTempTimes <- rep(NA, nreps)
  othTempLikes <- rep(NA, nreps)
  scale <- 0.0769827
  chol <- matrix(c(0.0488354,   -0.00871729,     1.52e-005, -3.48894e-005,
                   0,       0.49487, -2.15436e-005,  7.41979e-005,
                   0,             0,     0.0157383,  -0.000240194,
                   0,             0,             0,     0.0158859), nrow=4, ncol=4, byrow=T)
  
  propmat <- t(scale*chol)%*%(scale*chol)
  
  
  for(m in mvals){
    
    Rmodel <- nimbleModel(code, constants, data, inits)
    compileNimble(Rmodel)
    
    mcmcspec <- configureMCMC(Rmodel)
    mcmcspec$addSampler(type = 'RW_PFilter',
                        target = c('a', 'b', 'sigPN', 'sigOE'),
                        control = list(adaptive=F,
                                       propCov= propmat,
                                       adaptScaleOnly=F,
                                       m = m,
                                       latents = 'x'))
    
    
    mcmcspec$removeSamplers(ind = 1:104)
    Rmcmc <- buildMCMC(mcmcspec)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    for(r in rvals){
      for(n in 1:nreps){
        print(c(m,r,n))
        othTempTimes[n] <- sum(system.time(otherFilt(m,r))[1:2])

        
        nimTempTimes[n] <- sum(system.time(Cmcmc$run(r))[1:2])
      }    
      
      timelist <- rbind(timelist,cbind("Method" = c(paste("NIMBLE",filtType, sep=" "),
                                                    paste(otherPack,filtType,sep=" ")),
                                       "N. Particles" =c(m,m),
                                       "N. Iterations" = c(r,r),
                                        "Median" = c(median(nimTempTimes), median(othTempTimes))))
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
