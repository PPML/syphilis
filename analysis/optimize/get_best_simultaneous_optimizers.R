# library(syphLAMA)
devtools::load_all()
library(here)

# Find optimization.R output files
files_path <- here('inst/optims/7-19-19/')

# list the files
files <- grep("optim_list", list.files(files_path, full.names=T), value=T)

# read each file into a list
optims <- lapply(files, readRDS)


# Write helpers to extract data - posterior values and parameters
get_optim_vals <- function(optim_chain) {
  sapply(optim_chain, `[[`, 'value')
}

get_optim_pars <- function(optim_chain) {
  lapply(optim_chain, `[[`, 'par')
}

# extract posterior values
optim_vals <- lapply(optims, get_optim_vals)

# Find the best posterior values from each optim chain
maxs <- sapply(optim_vals, max) 
max_idxs <- sapply(optim_vals, function(x) which(x==max(x))[[1]])

# helper function to extract the last element of a list
last <- function(x) { x[[length(x)]] }

# Get the last optim output from each optim_chain
last_optims <- lapply(optims, last)
last_optims <- lapply(last_optims, function(x) x[['par']])
optim_pars <- do.call(rbind.data.frame, last_optims) # rbind to a dataframe

# Give optim_pars colnames
colnames(optim_pars) <- rep(names(theta_la), 2)

# Compute which parameters are the natural history parameters
natural_param_idxs <- which(colnames(optim_pars) %in% natural_history_parameters)

# Distinguish the first and second set of natural history parameters/columns of parameters.
first_natural_history_idxs <- natural_param_idxs[1:(length(natural_param_idxs)/2)]
second_natural_history_idxs <- natural_param_idxs[(length(natural_param_idxs)/2 + 1):length(natural_param_idxs)]

# Overwrite the first set of natural history parameters with the second set. 
optim_pars[,first_natural_history_idxs] <- optim_pars[,second_natural_history_idxs]

# Save the top25 simultaneous optimizers.
saveRDS(optim_pars, here("inst/optims/7-19-19/top25_best_simultaneous_params.rds"))

# Extract & Save the Top 25 Parameters for Each Location
optim_pars_LA <- optim_pars[, 1:(ncol(optim_pars)/2)]
optim_pars_MA <- optim_pars[, (ncol(optim_pars)/2) + 1:(ncol(optim_pars)/2)]

colnames(optim_pars_LA) <- names(theta_la)
colnames(optim_pars_MA) <- names(theta_ma)

saveRDS(optim_pars_LA, file.path(files_path, 'la_25_best_pars.rds'))
saveRDS(optim_pars_MA, file.path(files_path, 'ma_25_best_pars.rds'))


# Confirm that natural history parameters have been rewritten properly 
identical(optim_pars_LA[, natural_history_parameters], optim_pars_MA[, natural_history_parameters])


# Save top 5 parameters 
top5_maxs <- sort(maxs, decreasing=T)[1:5]
top5_idxs <- which(maxs %in% top5_maxs)

optim_top5_pars_LA <- optim_pars_LA[top5_idxs,]
optim_top5_pars_MA <- optim_pars_MA[top5_idxs,]

saveRDS(optim_top5_pars_LA, file.path(files_path, 'la_top5_best_pars.rds'))
saveRDS(optim_top5_pars_MA, file.path(files_path, 'ma_top5_best_pars.rds'))

# Get the best simultaneous optimizer
best_idx <- which(maxs == max(maxs))
theta_simultaneous <- unlist(optim_pars[best_idx,])
theta_la <- unlist(optim_pars_LA[best_idx,])
theta_ma <- unlist(optim_pars_MA[best_idx,])

# Write them into the package
usethis::use_data(theta_simultaneous, overwrite=TRUE)
usethis::use_data(theta_la, overwrite=TRUE)
usethis::use_data(theta_ma, overwrite=TRUE)
