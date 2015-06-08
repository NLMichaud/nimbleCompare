##  Contains code to run bootstrap and auxiliary particle filters.
##  For each kind of filter, we have a build function (buildBootF, buildAuxF),
##  and step functions.  There are two step functions for each filter: 
##  one which ends with NS (only saves latent state samples from
##  most recent time point), one which ends with S (saves latent state 
##  samples from all time points).  Also has a function (normMean) which returns the mean of a normally
##  distributed node.



auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), thresh_num=double()) 
    returnType(double())
)


# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.

bootFStepNS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num=double()) {
    returnType(double())
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(ess, double())
    declare(llEst, double(1,m))
    for(i in 1:m) {
      if(notFirst) {  
        copy(mv, model, nodes = 'xs', nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = 'x', row = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        llEst[i] <- wts[i]*mv['wts',i][1]
      }
      else{
        llEst[i] <- wts[i]/m
      }
    }
    # Normalize weights and calculate effective sample size 
    wts <- wts/sum(wts)
    ess <- 1/sum(wts^2) 
    
    # Determine whether to resample by weights or not
    if(ess < thresh_num){
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, nodes = 'x', nodesTo = 'xs', row = ids[i], rowTo = i)
        mv['wts',i][1] <<- 1/m
      }
    }
    else{
      for(i in 1:m){
        copy(mv, mv, "x", "xs", i, i)
        mv['wts',i][1] <<- wts[i]
      }
    }
    return(log(sum(llEst)))
  },  where = getLoadingNamespace()
)


bootFStepS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    # Get names of xs node for current and previous time point (used in copy)
    thisXSName <- paste("xs[,", t, "]", sep = "") 
    prevXSName <- paste("xs[,", t-1, "]", sep = "")
    thisXName <- paste("x[,", t, "]", sep = "")
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num=double()) {
    returnType(double())
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(ess, double())
    declare(llEst, double(1,m))
    for(i in 1:m) {
      if(notFirst) {  
        copy(mv, model, nodes = prevXSName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = thisXName, row = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        llEst[i] <- wts[i]*mv['wts',i][1,(t-1)]
      }
      else{
        llEst[i] <- wts[i]/m
      }
    }
    
    # Normalize weights and calculate effective sample size 
    wts <- wts/sum(wts)
    ess <- 1/sum(wts^2) 
    
    # Determine whether to resample by weights or not
    if(ess < thresh_num){
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, nodes = thisXName, nodesTo = thisXSName, row = ids[i], rowTo = i)
        mv['wts',i][1,t] <<- 1/m
      }
    }
    else{
      for(i in 1:m){
        copy(mv, mv, "x", "xs", i, i)
        mv['wts',i][1,t] <<- wts[i]
      }
    }
    return(log(sum(llEst)))
  },  where = getLoadingNamespace()
)

#' Creates a particle filter (sequential monte carlo) algorithm to estimate the log-likelihood for a model sub-graph
#'
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param thresh A number between 0 and 1 specifying when to resample: the resampling step will occur when the effective sample size is less than thresh*(number of particles)
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @author Daniel Turek
#' @details The resulting specialized particle filter algorthm will accept a single integer argument (m, default 10,000), which specifies the number of random \'particles\' to use for estimating the log-likelihood.  The algorithm returns the estimated log-likelihood value, and saves samples from the posterior distribution of the latent states in mv['xs',]
#' @examples
#' model <- nimbleModel(code = ...)
#' my_PF <- buildPF(model, 'x[1:100]')
#' Cmodel <- compileNimble(model)
#' Cmy_PF <- compileNimble(my_PF, project = model)
#' logLike <- Cmy_PF(m = 100000)
#' boot_X <- Cmy_PF$mv['xs',]

buildBootF <- nimbleFunction(
  setup = function(model, nodes, thresh = 0.5, silent = FALSE, saveAll = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
    if(0>thresh || 1<thresh || !is.numeric(thresh)) stop('thresh must be between 0 and 1')
    
    # Create mv variables for x state and sampled x states. 
    # Size depends on whethersaveAll=T or not 
    if(!saveAll){
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = dims[[1]],
                                                    xs = dims[[1]],
                                                    wts = dims[[1]])))
      bootStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes)){
        bootStepFunctions[[iNode]] <- bootFStepNS(model, mv, nodes, 
                                                  iNode, silent)  
      }
    }

    else{
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs', 'wts'),
                                        type = c('double', 'double', 'double'),
                                        size = list(x = c(dims[[1]],
                                                           length(dims)),
                                                    xs = c(dims[[1]],
                                                           length(dims)),
                                                    wts = c(dims[[1]],
                                                           length(dims))
                                                    )))
      bootStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes)){
        bootStepFunctions[[iNode]] <- bootFStepS(model, mv, nodes, 
                                                 iNode, silent)   
      }
    }
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    my_initializeModel$run()
    resize(mv, m)
    thresh_num <- ceiling(thresh*m)
    logL <- 0    
    for(iNode in seq_along(bootStepFunctions)) { 
      logL <- logL + bootStepFunctions[[iNode]]$run(m,thresh_num)
      if(logL == -Inf) return(logL)
    }
    return(logL)
  },  where = getLoadingNamespace()
)



auxFStepNS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    getmean <- normMean(model, thisNode)
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num=double()) {
    returnType(double())
    declare(auxWts, double(1,m))
    declare(auxl, double(1,m))
    declare(wts, double(1,m))
    declare(l, double(1, m))
    declare(ids, integer(1, m))
    ess <- 0
    if(notFirst){ #can't do auxiliary step for first time point
      for(i in 1:m) {
        copy(mv, model, 'x', prevNode, row=i)
        calculate(model, prevDeterm) 
        model[[thisNode]] <<- getmean$return_mean()  # returns E(x_t+1 | x_t)
        calculate(model, thisDeterm)
        auxl[i] <- exp(calculate(model, thisData))
        auxWts[i] <- auxl[i]*mv['wts',i][1]
      }
      auxWts <- auxWts/sum(auxWts)
      ess <- 1/sum(auxWts^2)
      if(ess<thresh_num){
        resamp <- 1
        rankSample(auxWts, m, ids, silent)
        for(i in 1:m){
          copy(mv, mv, 'x', 'xs', ids[i],i)
        }
      }
      else{
        resamp <- 0
        for(i in 1:m){
          copy(mv, mv, 'x', 'xs', i,i)
        }
      }
    }  
    for(i in 1:m) {
      if(notFirst) {  
        copy(mv, model, nodes = 'xs', nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = 'x', row = i)
      calculate(model, thisDeterm)
      l[i]  <- exp(calculate(model, thisData)) #likelihood
      #rescale weights by pre-sampling weight
      if(notFirst){
        if(resamp == 1){
          mv['wts',i][1] <<- l[i]/auxl[ids[i]]
        }
        else{
          mv['wts',i][1] <<- l[i]/auxl[i]
        }
      }
      else {
        mv['wts',i][1] <<- l[i]
      }
    }   
   
    # Skipping second resample step as per Carpenter 1999 / Cappe 06,
    # except for last time point.
    if(isLast){
      for(i in 1:m){
        wts[i] <- mv['wts',i][1]
      }
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, 'x', 'xs', ids[i],i)
      }
    }

    return(log(mean(l)))
  },  where = getLoadingNamespace()
)


auxFStepS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE) 
    #current time point
    t <- iNode
    # Get names of x and xs node for current and previous time point (used in copy)
    prevXSName <- paste("xs[,",t-1,"]",sep="")    
    thisXSName <- paste("xs[,",t,"]",sep="")
    prevXName <- paste("x[,",t-1,"]",sep="")
    thisXName <- paste("x[,",t,"]",sep="")
    getmean <- normMean(model, thisNode)
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num = double()) {
    returnType(double())
    declare(auxl, double(1,m))
    declare(auxWts, double(1,m))
    declare(wts, double(1,m))
    declare(ids, integer(1, m))
    declare(l, double(1,m))
    ess <- 0
    if(notFirst){
      for(i in 1:m) {
        copy(mv, model, prevXName, prevNode, row=i)
        calculate(model, prevDeterm) 
        model[[thisNode]] <<- getmean$return_mean() # returns E(x_t+1 | x_t)
        calculate(model, thisDeterm)
        auxl[i] <- exp(calculate(model, thisData))
        auxWts[i] <- auxl[i]*mv['wts',i][1,t-1]
      }
      auxWts <- auxWts/sum(auxWts)
      ess <- 1/sum(auxWts^2)
      if(ess<thresh_num){
        resamp <- 1
        rankSample(auxWts, m, ids, silent)
        for(i in 1:m){
          copy(mv, mv, prevXName, prevXSName, ids[i],i)
        }
      }
      else{
        resamp <- 0
        for(i in 1:m){
          copy(mv, mv, prevXName, prevXSName, i,i)
        }
      }
    }   
    for(i in 1:m) {
      if(notFirst) {
        copy(mv, model, nodes = prevXSName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = thisXName, row = i)
      calculate(model, thisDeterm)
      l[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        if(resamp == 1){
          mv['wts',i][1,t] <<- l[i]/auxl[ids[i]]
        }
        else{
          mv['wts',i][1,t] <<- l[i]/auxl[i]
        }
      }
      else{
        mv['wts',i][1,t] <<- l[i]
      }
    }
    
    if(isLast){
      for(i in 1:m){
        wts[i] <-  mv['wts',i][1,t] 
      }
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, thisXName, thisXSName, ids[i], i)
      }
    }
    return(log(mean(l)))
  },  where = getLoadingNamespace()
)


#' Creates an auxiliary particle filter 
#' 
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the Auxiliary particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param thresh A number between 0 and 1 specifying when to resample: the resampling step will occur when the effective sample size is less than thresh*(number of particles)
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @author Nick Michaud
#' @details Uses the Auxiliary Particle filter (Pitt & Shephard, 1999) algorithm which provides samples from f(x_t|y_t, x_t-1).  Sampled X values are stored in mv['xs',].  NOTE:  Filter currently assumes a normal transition equation f(x_t|x_t-1)!  Future functionality could include choice of g(x^(i)_t) function, currently restricted to E(x_t|x_t-1).
#' @examples
#' model <- nimbleModel(code = ...)
#' my_AuxF <- buildAuxF(model, 'x[1:100]', .9)
#' Cmodel <- compileNimble(model)
#' Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' logLike <- Cmy_AuxF(m = 100000)
buildAuxF <- nimbleFunction(
  setup = function(model, nodes, thresh=.5, silent = FALSE, saveAll = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) 
      stop('sizes or dimension of latent states varies')
    if(thresh<0 || thresh>1 || !is.numeric(thresh)) 
      stop('thresh must be between 0 and 1')
    
    #  Create mv variables for x state and sampled x states.  If saveLatent=T, the sampled x states
    #  will be recorded at each time point. 
    if(!saveAll){
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = dims[[1]],
                                                    xs = dims[[1]],
                                                    wts = dims[[1]])))
      auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes))
        auxStepFunctions[[iNode]] <- auxFStepNS(model, mv, nodes,
                                                iNode, silent) 
    }
    else{
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = c(dims[[1]],
                                                          length(dims)),
                                                    xs = c(dims[[1]],
                                                           length(dims)),
                                                    wts=c(dims[[1]],
                                                          length(dims)))))
      auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes))
        auxStepFunctions[[iNode]] <- auxFStepS(model, mv, nodes,
                                               iNode, silent) 
    }
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mv, m)  
    thresh_num <- ceiling(thresh*m)
    logL <- 0
    for(iNode in seq_along(auxStepFunctions)) {
      logL <- logL + auxStepFunctions[[iNode]]$run(m, thresh_num)      
      if(logL == -Inf) return(logL) }
    return(logL)
  },  where = getLoadingNamespace()
)

# Has a return_mean method which returns the mean 
# of a normally distributed nimble node.
normMean <- nimbleFunction(
  setup = function(model, node) {
    nfList <- nimbleFunctionList(node_stoch_dnorm)
    nfList[[1]] <- model$nodeFunctions[[node]]
  },  where = getLoadingNamespace(),
  methods = list(                        
    return_mean = function() {         
      returnType(double())           
      return(nfList[[1]]$get_mean()) 
    }                                  
  )                                      
)
