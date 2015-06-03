#m and thresh defined externally

#now let's code in SMC package
nStreams <- m
nPeriods <- length(data$y)
dimPerPeriod <- 1
threshnum <- m*thresh

#auxiliary proposals
generateStreamRepsFunc <- function(currentPeriod, lag1Streams, 
                                   lag1LogWeights, streamIndices ,
                                    ...){
  outMat <- matrix(nrow = m, ncol=1)
  if(currentPeriod==1){
    outMat[,1] <- lag1Streams[,1]
  }
  else{
    outMat[,1] <- inits$a*lag1Streams[,1]+inits$b
  }
  return(outMat)
}

#propogation step
generateNextStreamsFunc <- function(currentPeriod, lag1Streams, 
                                    lag1LogWeights, streamIndices, 
                                    streamReps, startingStreams, ...){
  outMat <- matrix(nrow = m, ncol=1)
  if(currentPeriod==1){
    outMat[,1]<- rnorm(m, inits$b/(1-inits$a),  sqrt(inits$sigPN^2/(1-inits$a^2)))
  }
  else{
    outMat[,1] <- rnorm(m, lag1Streams[streamIndices,1]*inits$a + inits$b, inits$sigPN)
  }
  return(outMat)
}


#calculate f(y_t|x_t)
logObsDensFunc <- function(currentPeriod, currentStreams, ...){
  return(dnorm(data$y[currentPeriod], mean = currentStreams[,1], sd=inits$sigOE,log=TRUE))
}

#always resample for testing purposes
resampCriterionFunc <- function(currentPeriod, currentStreams, currentLogWeights, ...){
  #wts <- exp(currentLogWeights)
  #wts <- wts/sum(wts)
  #ess <- 1/sum(wts^2)
  #if(ess > threshnum) return(FALSE)
  #else return(TRUE)
  return(TRUE)  
}

smcauxfunc <- function(nps){auxiliaryParticleFilter(nps, nPeriods, dimPerPeriod,
                                                  generateStreamRepsFunc, generateNextStreamsFunc,
                                                  logObsDensFunc, resampCriterionFunc, 
                                                  resampFunc=NULL, summaryFunc=NULL, nMHSteps=0,
                                                  returnStreams       = TRUE, 
                                                  returnLogWeights    = TRUE, 
                                                  startingStreams=NULL,
                                                  nStreamsPreResamp=nps,
                                                  verboseLevel        = 0,     
                                                  inits, data, threshnum,nps)}

