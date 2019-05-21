
# In this document, we will continue the optimization of the dLogPosterior 
# function.

# dependencies
library(here)
devtools::load_all()

# config
N_loops <- 20
input_directory <- here("inst/optims/4-8-19")
output_directory <- here("inst/optims/5-21-19")

# command line arguments
state <<- commandArgs(trailingOnly=T)[[1]]
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[2]])

# load model
load.start()

# read old optim chain
optim_chain <- 
    readRDS(file.path(
              input_directory,
              paste0(state,
              "_optim_list_",
              nth_sim,
              ".rds")
            ))

# get the original theta vector names properly 
orig_theta_names <- names(get(paste0("theta_", tolower(state))))

# write a version of dLogPosterior that names the theta vector
# since the optim() function calls its argument f with 
# an unnamed numeric vector.
dLogPosterior_for_optim <- function(theta) {
	names(theta) <- orig_theta_names
	d <- dLogPosterior(theta)
	print(d)
	return(d)
}


# Write helpers to extract prior optimization data
get_optim_vals <- function(optim_chain) {
  sapply(optim_chain, `[[`, 'value')
}

# Extract posterior values
optim_vals <- sapply(1:length(optim_chain), function(x) dLogPosterior_for_optim(optim_chain[[x]]$par))
max_idx <- which(optim_vals == max(optim_vals))
theta <- orig_theta <- optim_chain[[max_idx]]$par

# dLogPosterior_val <- dLogPosterior_for_optim(theta)

print(paste0("starting with dLogPosterior_val: ", max(optim_vals)))
print(paste0("confirming dLogPosterior_val matches:", dLogPosterior_val))


# main optim loop -- run optims until 'done' and checkpoint outputs iteratively 
n_done <- 0
done <- F
optim_list <- list()
while (!done) {
  out <- optim(
    par = theta,
    fn = dLogPosterior_for_optim,
    method = c("Nelder-Mead", "BFGS")[[if (n_done > 5) 1 else 2]],
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
              site,
              "_optim_list_",
              nth_sim,
              ".rds"
            )
    )
  } else {
    theta <- orig_theta
 }
}
