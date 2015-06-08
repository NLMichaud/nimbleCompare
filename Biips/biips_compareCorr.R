
biipsDat <- constants
biipsDat$y <- data$y
biipsDat$a <- inits$a
biipsDat$b <- inits$b
biipsDat$sigPN <- inits$sigPN
biipsDat$sigOE <- inits$sigOE
bmod <- biips_model("corrMod.BUGS", biipsDat, quiet=TRUE)
bipLikeFunc  <- function(m){ biips_smc_samples(bmod, 'x[1:100]', m, type = "f",
                                               rs_thres=1 )$log_marg_like}