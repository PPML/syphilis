
# In this document, we will continue the optimization of the dLogPosterior 
# function.

# dependencies
library(here)
devtools::load_all()

# config
N_loops <- 20
input_directory <- here("inst/optims/10-11-19/")
output_directory <- here("inst/optims/10-12-19/")

# command line arguments
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[1]])

### Construct dLogPosterior_simultaneous to calibrate to both 
source(here("analysis/optimize/construct_simultaneous_optimization_function.R"))

# Read 
top_optim_trace_path_la <- paste0(input_directory, "la_25_best_pars.rds")
optim_trace_la <- readRDS(top_optim_trace_path_la)
theta_la <- unlist(optim_trace_la[as.numeric(nth_sim), ])

top_optim_trace_path_ma <- paste0(input_directory, "ma_25_best_pars.rds")
optim_trace_ma <- readRDS(top_optim_trace_path_ma)
theta_ma <- unlist(optim_trace_ma[as.numeric(nth_sim), ])

theta <- c(theta_la, theta_ma)
orig_theta <- theta 

dLogPosterior_val <- dLogPosterior_simultaneous(theta)

print(paste0("starting with dLogPosterior_simultaneous value: ", dLogPosterior_val))

# main optim loop -- run optims until 'done' and checkpoint outputs iteratively 
n_done <- 0
done <- F
optim_list <- list()
while (n_done <= N_loops) {
  out <- optim(
    par = theta,
    fn = dLogPosterior_simultaneous,
    method = c("Nelder-Mead", "BFGS")[[if (n_done > 10) 1 else 2]],
    control = list(fnscale = -1, maxit=1000)
  )

  if (out$value != -1e32) {
		# update and save last optim loop
    n_done <- n_done + 1
    done = n_done == N_loops
    optim_list[[n_done]] <- out
		out$value
		if (out$value > dLogPosterior_val) {
			theta <- out$par
			dLogPosterior_val <- out$value
		}
    saveRDS(optim_list,
            paste0(
              output_directory,
              "optim_list_",
              nth_sim,
              ".rds"
            )
    )
  } else {
    theta <- orig_theta
 }
}
