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

