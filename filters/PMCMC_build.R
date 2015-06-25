###################################################################################
### RW_llFunction, does a block RW, but using a generic log-likelihood function ###
###################################################################################

sampler_RW_llFunctionBlock <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ###  control list extraction  ###
    adaptive       <- control$adaptive
    adaptScaleOnly <- control$adaptScaleOnly
    adaptInterval  <- control$adaptInterval
    scale          <- control$scale
    propCov        <- control$propCov 
    llFunction     <- control$llFunction
    includesTarget <- control$includesTarget
  
    ###  node list generation  ###
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    ###  numeric value generation  ###
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory          <- c(0, 0)
    acceptanceRateHistory <- c(0, 0)
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
    if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    statSums  <- matrix(0, nrow=1, ncol=d)   # sums of each node, stored as a row-matrix
    statProds <- matrix(0, nrow=d, ncol=d)   # sums of pairwise products of nodes
    ###  nested function and function list definitions  ###
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d)
  },
  
  run = function() {
    modelLP0 <- llFunction$run()
    if(!includesTarget)     modelLP0 <- modelLP0 + getLogProb(model, target)
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    modelLP1 <- llFunction$run()+ getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(adaptive)     adaptiveProcedure(jump)
  },
  
  methods = list(
    
    generateProposalVector = function() {
      ##declare(origValueVector, double(1, d))
      ##origValueVector <- values(model, target)
      ##declare(normalVarVector, double(1, d))
      ##for(i in 1:d)     {   normalVarVector[i] <- rnorm(1, 0, 1)   }
      ##propValueMatrix <- asRow(origValueVector) + asRow(normalVarVector) %*% chol_propCov * scale
      ##propValueVector <- propValueMatrix[1, ]
      propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov * scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly) {
        declare(newValues, double(1, d))
        
        newValues <- values(model, target)
        statSums  <<- statSums + asRow(newValues)
        statProds <<- statProds + asCol(newValues) %*% asRow(newValues)
      }
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        setSize(scaleHistory,          timesAdapted)
        setSize(acceptanceRateHistory, timesAdapted)
        scaleHistory[timesAdapted] <<- scale
        acceptanceRateHistory[timesAdapted] <<- acceptanceRate
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$gamma1
          empirCov <- (statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
          statSums  <<- statSums  * 0
          statProds <<- statProds * 0      ##  setAll(statProds, 0)    ## setAll() doesn't work in R, and doesn't work for vectors (only works for dim=2 objects)
        }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      scaleHistory          <<- scaleHistory          * 0
      acceptanceRateHistory <<- acceptanceRateHistory * 0
      statSums  <<- statSums  * 0
      statProds <<- statProds * 0
      my_calcAdaptationFactor$reset()
    }
  ),  where = getLoadingNamespace()
)

#######################################################################################  
### RW_llFunction, does a block RW, but using a particle filter likelihood function ###
#######################################################################################

sampler_RW_PFilter <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ###  control list extraction  ###
    adaptive       <- control$adaptive
    adaptScaleOnly <- control$adaptScaleOnly
    adaptInterval  <- control$adaptInterval
    scale          <- control$scale
    propCov        <- control$propCov 
    m              <- control$m
    latents        <- control$latents
    if(is.na(m)){
      m <- 1000
    }
    
    
    
    ###  node list generation  ###
    
    if(!all(target%in%Rmodel$getNodeNames(stochOnly=T, includeData=F,
                                     topOnly=T))){
      stop("PMCMC target can only consist of top level parameters")
    }
    
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    
    calcNodes <- model$getDependencies(target)
    ###  numeric value generation  ###
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    scaleHistory          <- c(0, 0)
    acceptanceRateHistory <- c(0, 0)
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
    if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    #if(!is.integer(m))            stop('m must be an integer')
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    statSums  <- matrix(0, nrow=1, ncol=d)   # sums of each node, stored as a row-matrix
    statProds <- matrix(0, nrow=d, ncol=d)   # sums of pairwise products of nodes

    storeLP <- modelValues(modelValuesSpec(vars = c('LP0'),
                                      type = c('double'),
                                      size = list(LP0 = 1)
                                      ))
    ###  nested function and function list definitions  ###
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d)
    my_particleFilter <- buildBootF(model, latents, .5, saveAll=F)
  },
  
  run = function() {
    if(is.nan(storeLP['LP0',1][1])){
      modelLP0 <- my_particleFilter$run(m)
      storeLP['LP0',1][1] <<- modelLP0 + getLogProb(model, target)
    }
    modelLP0 <- storeLP['LP0',1][1]
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    modelLP1 <- my_particleFilter$run(m)
    modelLP1 <- modelLP1 + getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(jump)         storeLP['LP0',1][1] <<- modelLP1
    if(adaptive)     adaptiveProcedure(jump)
  },
  
  methods = list(
    
    generateProposalVector = function() {
      ##declare(origValueVector, double(1, d))
      ##origValueVector <- values(model, target)
      ##declare(normalVarVector, double(1, d))
      ##for(i in 1:d)     {   normalVarVector[i] <- rnorm(1, 0, 1)   }
      ##propValueMatrix <- asRow(origValueVector) + asRow(normalVarVector) %*% chol_propCov * scale
      ##propValueVector <- propValueMatrix[1, ]
      propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov * scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly) {
        declare(newValues, double(1, d))
        
        newValues <- values(model, target)
        statSums  <<- statSums + asRow(newValues)
        statProds <<- statProds + asCol(newValues) %*% asRow(newValues)
      }
      if(timesRan %% adaptInterval == 0) {

        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        setSize(scaleHistory,          timesAdapted)
        setSize(acceptanceRateHistory, timesAdapted)
        scaleHistory[timesAdapted] <<- scale
        acceptanceRateHistory[timesAdapted] <<- acceptanceRate
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$gamma1
          empirCov <- (statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
          statSums  <<- statSums  * 0
          statProds <<- statProds * 0      ##  setAll(statProds, 0)    ## setAll() doesn't work in R, and doesn't work for vectors (only works for dim=2 objects)
        }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      scaleHistory          <<- scaleHistory          * 0
      acceptanceRateHistory <<- acceptanceRateHistory * 0
      statSums  <<- statSums  * 0
      statProds <<- statProds * 0
      my_calcAdaptationFactor$reset()
    }
  ),  where = getLoadingNamespace()
)


sampler_RW_PFilterResamp <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ###  control list extraction  ###
    adaptive       <- control$adaptive
    adaptScaleOnly <- control$adaptScaleOnly
    adaptInterval  <- control$adaptInterval
    scale          <- control$scale
    propCov        <- control$propCov 
    m              <- control$m
    latents        <- control$latents
    if(is.na(m)){
      m <- 1000
    }
    
    ###  node list generation  ###
    
    if(!all(target%in%model$getNodeNames(stochOnly=T, includeData=F,
                                         topOnly=T))){
      stop("PMCMC target can only consist of top level parameters")
    }
    
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    
    calcNodes <- model$getDependencies(target)
    ###  numeric value generation  ###
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory          <- c(0, 0)
    acceptanceRateHistory <- c(0, 0)
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
    if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    #if(!is.integer(m))            stop('m must be an integer')
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    statSums  <- matrix(0, nrow=1, ncol=d)   # sums of each node, stored as a row-matrix
    statProds <- matrix(0, nrow=d, ncol=d)   # sums of pairwise products of nodes
    ###  nested function and function list definitions  ###
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d)
    my_particleFilter <- buildBootF(model, latents, .5, , saveAll=F)
  },
  
  run = function() {
    modelLP0 <- my_particleFilter$run(m)
    modelLP0 <- modelLP0 + getLogProb(model, target)
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    modelLP1 <- my_particleFilter$run(m)
    modelLP1 <- modelLP1 + getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(adaptive)     adaptiveProcedure(jump)
  },
  
  methods = list(
    
    generateProposalVector = function() {
      ##declare(origValueVector, double(1, d))
      ##origValueVector <- values(model, target)
      ##declare(normalVarVector, double(1, d))
      ##for(i in 1:d)     {   normalVarVector[i] <- rnorm(1, 0, 1)   }
      ##propValueMatrix <- asRow(origValueVector) + asRow(normalVarVector) %*% chol_propCov * scale
      ##propValueVector <- propValueMatrix[1, ]
      propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov * scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly) {
        declare(newValues, double(1, d))
        
        newValues <- values(model, target)
        statSums  <<- statSums + asRow(newValues)
        statProds <<- statProds + asCol(newValues) %*% asRow(newValues)
      }
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        setSize(scaleHistory,          timesAdapted)
        setSize(acceptanceRateHistory, timesAdapted)
        scaleHistory[timesAdapted] <<- scale
        acceptanceRateHistory[timesAdapted] <<- acceptanceRate
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$gamma1
          empirCov <- (statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
          statSums  <<- statSums  * 0
          statProds <<- statProds * 0      ##  setAll(statProds, 0)    ## setAll() doesn't work in R, and doesn't work for vectors (only works for dim=2 objects)
        }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      scaleHistory          <<- scaleHistory          * 0
      acceptanceRateHistory <<- acceptanceRateHistory * 0
      statSums  <<- statSums  * 0
      statProds <<- statProds * 0
      my_calcAdaptationFactor$reset()
    }
  ),  where = getLoadingNamespace()
)

