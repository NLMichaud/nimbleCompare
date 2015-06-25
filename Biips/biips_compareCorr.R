
biipsDat <- constants
biipsDat$y <- data$y
biipsDat$a <- inits$a
biipsDat$b <- inits$b
biipsDat$sigPN <- inits$sigPN
biipsDat$sigOE <- inits$sigOE
bmod <- biips_model("corrMod.BUGS", biipsDat, quiet=TRUE)
bipLikeFunc  <- function(m){ biips_smc_samples(bmod, 'x[1:100]', m, type = "f",
                                               rs_thres=1 )$log_marg_like}
bipPMCMC <- function(m,r){
  bpmcmc <- biips_pmmh_init(bmod, param_names = c("a","b","sigOE","sigPN"), latent_names = "x[1:100]",
                            inits =  list(a = 0.95, b = 1, sigOE = 0.05, sigPN = 0.2),
                            rw_step = list(a = 0.000014, b = .00145, sigOE = 0.0000014, sigPN = 0.0000015),
                            transform=F,
                            n_rescale = 200)
  return(biips_pmmh_samples(bpmcmc, r, m,rs_thres=1 ))
}


