##Produce tables and plots for PMCMC comparison

require(nimble)
require(ggplot2)
require(mvtnorm)
require(pomp, quietly=TRUE)
require(coda)
require(Rbiips)

load("data/model_SSMcorrelated.RData")
source("filters/bootF_build.R")
source("filters/PMCMC_build.R")


thresh <- 1

#define model in different packages
source("pomp/pomp_compareCorr.R", echo=F, verbose=F)
source("Biips/biips_compareCorr.R", echo=F, verbose=F)

# Helpful functions for making comparison tables and plots
source("compare_help_functions.R") 

scale <- 0.0769827
chol <- matrix(c(0.0488354,   -0.00871729,     1.52e-005, -3.48894e-005,
                 0,       0.49487, -2.15436e-005,  7.41979e-005,
                 0,             0,     0.0157383,  -0.000240194,
                 0,             0,             0,     0.0158859), nrow=4, ncol=4, byrow=T)

propmat <- t(scale*chol)%*%(scale*chol)

nMeths <- 6


methNms <- c("NIMBLE PMCMC Resamp Adapt","NIMBLE PMCMC Resamp NoAdapt","NIMBLE PMCMC NoResamp Adapt",
             "NIMBLE PMCMC NoResamp NoAdapt", "POMP PMCMC", "Biips PMCMC")
methFuncs <- function(nm, m){
  if(nm == "NIMBLE PMCMC Resamp Adapt"){
    Rmodel <- nimbleModel(code, constants, data, inits)
    compileNimble(Rmodel)
    
    mcmcspec <- configureMCMC(Rmodel)
    mcmcspec$addSampler(type = 'RW_PFilterResamp',
                        target = c('a', 'b', 'sigPN', 'sigOE'),
                        control = list(adaptive=T,
                                       propCov= propmat,
                                       adaptScaleOnly=F,
                                       m = m,
                                       latents = 'x'))
    
    
    mcmcspec$removeSamplers(ind = 1:104)
    Rmcmc <- buildMCMC(mcmcspec)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    return(Cmcmc)
  }
  
  if(nm == "NIMBLE PMCMC Resamp NoAdapt"){
    Rmodel <- nimbleModel(code, constants, data, inits)
    compileNimble(Rmodel)
    
    mcmcspec <- configureMCMC(Rmodel)
    mcmcspec$addSampler(type = 'RW_PFilterResamp',
                        target = c('a', 'b', 'sigPN', 'sigOE'),
                        control = list(adaptive=F,
                                       propCov= propmat,
                                       adaptScaleOnly=F,
                                       m = m,
                                       latents = 'x'))
    
    
    mcmcspec$removeSamplers(ind = 1:104)
    Rmcmc <- buildMCMC(mcmcspec)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    return(Cmcmc)
  }
  
  if(nm == "NIMBLE PMCMC NoResamp Adapt"){
    Rmodel <- nimbleModel(code, constants, data, inits)
    compileNimble(Rmodel)
    
    mcmcspec <- configureMCMC(Rmodel)
    mcmcspec$addSampler(type = 'RW_PFilter',
                        target = c('a', 'b', 'sigPN', 'sigOE'),
                        control = list(adaptive=T,
                                       propCov= propmat,
                                       adaptScaleOnly=F,
                                       m = m,
                                       latents = 'x'))
    
    
    mcmcspec$removeSamplers(ind = 1:104)
    Rmcmc <- buildMCMC(mcmcspec)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    return(Cmcmc)
  }
  
  if(nm == "NIMBLE PMCMC NoResamp NoAdapt"){
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
    return(Cmcmc)
  }
  
  if(nm == "POMP PMCMC"){
    return(pompMCMC)
  } 
  if(nm == "Biips PMCMC"){
    return(bipPMCMC)
  }
}

timelist <- data.frame("Method" = NA,
                       "N. Particles" =NA,
                       "N. Iterations "= NA,
                       "Median" = NA)

nreps = 5
tempTimes <- rep(NA, nreps)



mvals <- c(100, 500, 1000)
rvals <- c(100, 1000) 
for(i in 1:nMeths){
  methNm <- methNms[i]
  for(m in mvals){
    meth <- methFuncs(methNm, m)    
    for(r in rvals){
      for(n in 1:nreps){ 
        if(!(methNm %in% c("POMP PMCMC", "Biips PMCMC"))){
          tempTimes[n] <- sum(system.time(meth$run(r, reset = TRUE))[1:2])
        }
        else{
          tempTimes[n] <- sum(system.time( meth(m,r))[1:2])          
        }
      }
      gc()
      timelist <- rbind(timelist, c("Method" = methNm,
                                    "N. Particles" =m,
                                    "N. Iterations" = r,
                                    "Median" = median(tempTimes)))
      save(timelist, file="timeList.RData")
    }
    gc()
  }
}


burnIn <- 2500

m <- 1000
r <- 7500

essTable <- data.frame("Method" = NA, "a" = NA, "b" = NA, "sigOE" = NA, "sigPN" = NA)

#run for all packages except Biips
for(i in 1:(length(methNms)-1)){
  methNm <- methNms[i]
  meth <- methFuncs(methNm, m)
  
  if(!(methNm %in% c("POMP PMCMC", "Biips PMCMC"))){
    meth$run(r)
    tempDat <- as.mcmc(as.matrix(meth$mvSamples)[-(1:burnIn),])
  }
  else{
    tempDat <- meth(m,r)
    tempDat <- as.mcmc(conv.rec(tempDat)[-(1:burnIn),4:7])
  }
  
  essTable <- rbind(essTable, c("Method" = methNm, effectiveSize(tempDat) ))
  save(essTable, file="essTable.Rdata")
  png(paste("plots//pjlot", i, ".png",sep=""))
  plot(tempDat, ask=F)
  dev.off()
}

