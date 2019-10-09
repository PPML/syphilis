library(here)
devtools::load_all()

# command line arguments
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[1]])
if (is.null(nth_sim)) nth_sim <- 1

# script config
N_loops <- 20
output_directory <- here("inst/optims/10-8-19/")

source(here("analysis/optimize/construct_simultaneous_optimization_function.R"))

# main optim loop -- run optims until 'done' and checkpoint outputs iteratively 
n_done <- 0
done <- F
optim_list <- list()
current_val <- dLogPosterior_simultaneous(theta)

while (!done) {
  out <- optim(
    par = theta,
    fn = dLogPosterior_simultaneous,
    method = "Nelder-Mead",
    control = list(fnscale = -1, maxit=1000)
  )

  if (out$value > -1e32 && out$value > current_val) {
		# update and save last optim loop
    n_done <- n_done + 1
    done = n_done == N_loops
    optim_list[[n_done]] <- out
    theta <- out$par
		current_val <- out$value
    saveRDS(optim_list,
            paste0(
              output_directory,
              "optim_list_",
              nth_sim,
              ".rds"
            )
    )
  }
}

