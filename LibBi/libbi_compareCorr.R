require("bi")
init_bi <- bi_wrapper$new(client="filter", model_file_name = "corrMod.bi",
                          config="model_config.conf")


bibootLike <- function(m){
  init_bi$run(add_options=paste("--init-file initvars.nc  --filter bootstrap  --nparticles ", m, sep=""))
  return(bi_read_var(init_bi$result$output_file_name, "LL"))
}

biauxLike <- function(m){
  init_bi$run(add_options=paste("--init-file initvars.nc  --filter lookahead  --nparticles ", m, sep=""))
  return(bi_read_var(init_bi$result$output_file_name, "LL"))
}