library(here)
devtools::load_all()

# command line arguments
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[1]])
if (is.null(nth_sim)) nth_sim <- 1

# script config
N_loops <- 20
output_directory <- here("inst/optims/6-6-19/")

# save orig_theta_names for use in dLogPosterior_by_state
orig_theta_names <- names(theta_la)

# save all of the global variable names loaded on package load
globals <- ls()

# construct a new environment and save everything into it for each state
for (state in c('LA', 'MA')) {
	# parametrize the model 
	load.start()

	# construct an environment for the state
	state_env <- new.env()

	# save each of the parametrizing variables into the new environment
	for (obj in setdiff(ls(), c(globals, 'nth_sim', 'orig_theta_names', 'state',
	'state_env', 'globals', 'la_env', 'ma_env', 'theta_la', 'theta_ma'))) { 
		assign(x = obj, value = get(obj), envir = state_env)
	}

	# export the state env with la_env or ma_env naming
	assign(x = paste0(tolower(state), '_env'), value = state_env)


	# get the best parametrization of that state
	# optim_trace <- readRDS(system.file('optims', paste0(tolower(state), '_top5_pars.rds'), package='syphLAMA'))
	optim_trace <- readRDS(here(paste0('inst/optims/', tolower(state), '_top5_pars.rds')))
	dLogPosterior_vals <- apply(optim_trace, 1, dLogPosterior)
	max_idx <- which(dLogPosterior_vals == max(dLogPosterior_vals))
	assign(x = paste0('theta_', tolower(state)), value = unlist(optim_trace[max_idx, ]))

	# remove all of the parametrizing variables except as they're stored in the state_env
	rm(list=setdiff(ls(), c(globals, 'orig_theta_names', 'globals', 'la_env', 'ma_env', 'theta_la', 'theta_ma')))
}


dLogPosterior_by_state <- function(theta, state) {
	state_env <- get(paste0(tolower(state), '_env'))
	attach(state_env) 
	names(theta) <- orig_theta_names
	val <- dLogPosterior(theta)
	detach(state_env)
	return(val)
}

theta <- c(theta_la, theta_ma)

dLogPosterior_simultaneous <- function(theta) {
	# overwrite the natural history parameters to be the same for each state in theta
	theta[1:(length(theta)/2)][natural_history_parameters] <- 
		theta[length(theta)/2 + 1:(length(theta)/2)][natural_history_parameters]

	la_val <- suppressMessages(dLogPosterior_by_state(theta[1:(length(theta)/2)], 'LA'))
	ma_val <- suppressMessages(dLogPosterior_by_state(theta[length(theta)/2 + 1:(length(theta)/2)], 'MA'))

	if (any(-1e20 >= c(la_val, ma_val)) { 
		# try again, for some reason dLogPosterior occasionally returns -1e20 for the same values
		la_val <- suppressMessages(dLogPosterior_by_state(theta[1:(length(theta)/2)], 'LA'))
		ma_val <- suppressMessages(dLogPosterior_by_state(theta[length(theta)/2 + 1:(length(theta)/2)], 'MA'))
	}
	return(la_val + ma_val)
}


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
              site,
              "_optim_list_",
              nth_sim,
              ".rds"
            )
    )
  }
}


dLogPosterior_simultaneous_print <- function(x) {
	val <- dLogPosterior_simultaneous(x)
	print(val)
	return(val)
}

  out <- optim(
    par = theta,
    fn = dLogPosterior_simultaneous_print,
    method = c("BFGS", "Nelder-Mead")[[if (n_done > 15) 1 else 2]],
    control = list(fnscale = -1, maxit=10)
  )
