# Update the theta_la and theta_ma vectors

devtools::load_all(".")
library(here)

# Load Data
# trace_la <- readRDS(here('inst/optims/6-6-19/la_top5_best_pars.rds'))
# trace_ma <- readRDS(here('inst/optims/6-6-19/ma_top5_best_pars.rds'))

# Compute the dLogPosterior values for Each Entry of Each trace row
# for (state in c('LA', 'MA')) { 
# 	load.start()
# 	assign(tolower(paste0(state, '_vals')),
# 		apply(get(paste0('trace_', tolower(state))), 1, dLogPosterior))
# }

# la_best_idx <- which(la_vals == max(la_vals))[[1]]
# ma_best_idx <- which(ma_vals == max(ma_vals))[[1]]

# theta_la <- unlist(trace_la[la_best_idx,])
# theta_ma <- unlist(trace_ma[ma_best_idx,])

# usethis::use_data(theta_la, overwrite = TRUE)
# usethis::use_data(theta_ma, overwrite = TRUE)

source(here("analysis/optimize/construct_simultaneous_optimization_function.R"))

optim_pars <- readRDS(here("inst/optims/6-6-19/top25_best_simultaneous_params.rds"))
post_vals <- apply(optim_pars, 1, 
	function(x) dLogPosterior_simultaneous(unlist(x)))


best_idx <- which(post_vals == max(post_vals))[[1]]

theta <- unlist(optim_pars[best_idx,])
names(theta) <- rep(names(theta_la), 2)

theta_la <- theta[1:(length(theta)/2)]
theta_ma <- theta[(length(theta)/2+1):length(theta)]

usethis::use_data(theta_la, overwrite = TRUE)
usethis::use_data(theta_ma, overwrite = TRUE)
