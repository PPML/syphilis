# Run AMCMC With Top Values from Optimization
library(adaptMCMC)
library(here)
devtools::load_all()

args <- commandArgs(trailingOnly = T)
array_id <- args[[1]]
output_dir <- here("inst/mcmc/7-1-19/")


### Construct dLogPosterior_simultaneous to calibrate to both 
### Louisiana and Massachusetts at the same time sharing the 
### natural_history_parameters.

# save orig_theta_names for use in dLogPosterior_by_state
data("theta_la")
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


dLogPosterior_simultaneous <- function(theta) {
	# overwrite the natural history parameters to be the same for each state in theta
	theta[1:(length(theta)/2)][natural_history_parameters] <- 
		theta[length(theta)/2 + 1:(length(theta)/2)][natural_history_parameters]

	la_val <- suppressMessages(dLogPosterior_by_state(theta[1:(length(theta)/2)], 'LA'))
	ma_val <- suppressMessages(dLogPosterior_by_state(theta[length(theta)/2 + 1:(length(theta)/2)], 'MA'))

	if (any(-1e20 >= c(la_val, ma_val))) { 
		# try again, for some reason dLogPosterior occasionally returns -1e20 for the same values
		la_val <- suppressMessages(dLogPosterior_by_state(theta[1:(length(theta)/2)], 'LA'))
		ma_val <- suppressMessages(dLogPosterior_by_state(theta[length(theta)/2 + 1:(length(theta)/2)], 'MA'))
	}
	return(la_val + ma_val)
}

top_optim_trace_path_la <- paste0(here(paste0("inst/optims/6-6-19/la_25_best_pars.rds")))
optim_trace_la <- readRDS(top_optim_trace_path_la)
theta_la <- unlist(optim_trace_la[ (as.numeric(array_id) %% 25)+1, ])

top_optim_trace_path_ma <- paste0(here(paste0("inst/optims/6-6-19/ma_25_best_pars.rds")))
optim_trace_ma <- readRDS(top_optim_trace_path_ma)
theta_ma <- unlist(optim_trace_ma[ (as.numeric(array_id) %% 25)+1, ])

sd_optims_la <- apply(optim_trace_la, 2, sd)
sd_optims_ma <- apply(optim_trace_ma, 2, sd)
sd_optims <- c(sd_optims_la, sd_optims_ma)

theta <- c(theta_la, theta_ma)

sample <- MCMC(
  p = dLogPosterior_simultaneous,
  n = 25000,
  init = theta,
  scale = sd_optims / 1000,
  adapt = T,
  acc.rate = .234, showProgressBar = T)

saveRDS(sample, file = paste0(output_dir, "mcmc_25000_", array_id, ".rds"))
