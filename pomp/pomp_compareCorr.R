mod.step<- Csnippet("
  double eps = rnorm(0,sigPN);
  x = a*x+b+eps;
")

rmeas <- Csnippet("
  y = rnorm(x,sigOE);
")

dmeas <- Csnippet("

    lik = dnorm(y,x,sigOE,give_log);
")

                   
initializer <- function(params,t0,...){
  x0 <- rnorm(1, params['b']/(1-params['a']), sqrt(as.numeric(params['sigPN'])^2 / (1-params['a']^2)))
  return(c(x=x0))
}

pompdat <- data.frame(y=data$y, times=c(1:100))
modelPomp <- pomp(data=pompdat, time="times", rprocess=discrete.time.sim(mod.step,delta.t=1),
                  dmeasure=dmeas, 
                  statenames="x",paramnames=c("b","a","sigPN", "sigOE"),
                  t0=1, initializer=initializer, obsnames="y", give_log=1)

pompPf <- function(nps){
  tmp <- pfilter(modelPomp,Np=nps,params=c(b=inits$b,a=inits$a,sigPN=inits$sigPN,sigOE=inits$sigOE),
              save.states=T)
  return(logLik(tmp))
}


#cant get rprior for pomp to work in C, so just simulate initial values in R :(
pompLW <- function(nps){
  simPrior <- matrix(nrow=4, ncol=nps)
  simPrior[1,] <- rnorm(nps, 0, sd = 1000)
  simPrior[2,] <- runif(nps, -.9999, .9999)
  simPrior[3,] <- runif(nps, 0.0001, 1)
  simPrior[4,] <- runif(nps, 0.0001, 1)
  rownames(simPrior) <- c("b","a","sigPN", "sigOE")
  tmp <- bsmc(modelPomp, params=simPrior, nP=nps)
  return(tmp$log.evidence)
}

pompMCMC <- function(nps, nreps){
  scale <- 0.0769827
  chol <- matrix(c(0.0488354,   -0.00871729,     1.52e-005, -3.48894e-005,
                   0,       0.49487, -2.15436e-005,  7.41979e-005,
                   0,             0,     0.0157383,  -0.000240194,
                   0,             0,             0,     0.0158859), nrow=4, ncol=4, byrow=T)
  
  propmat <- t(scale*chol)%*%(scale*chol)
  prop <- function (theta) 
  {
    thetaTemp <- rmvnorm(1, theta, propmat)
    
    while((thetaTemp[1]>1)|| (thetaTemp[3]<0) || (thetaTemp[4] < 0))   thetaTemp <- rmvnorm(1, theta, propmat)
    theta <- c("a" = thetaTemp[1], "b" = thetaTemp[2], "sigOE" = thetaTemp[3], "sigPN" = thetaTemp[4])
    theta
  }
  
  svals <- c(a = 0.95, b = 1, sigOE = 0.05, sigPN = 0.2)
 return( pmcmc(object = modelPomp, Nmcmc = nreps, start = svals, proposal = prop, Np = nps,
        max.fail = 10000000))
}

