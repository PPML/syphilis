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

# print(sort(maxs, decreasing=T))

last <- function(x) { x[[length(x)]] }

last_optims <- lapply(optims, last)
last_optims <- lapply(last_optims, function(x) x[['par']])
optim_pars <- do.call(rbind.data.frame, last_optims)

colnames(optim_pars) <- rep(names(theta_la), 2)
natural_param_idxs <- which(colnames(optim_pars) %in% natural_history_parameters)


first_natural_history_idxs <- natural_param_idxs[1:(length(natural_param_idxs)/2)]
second_natural_history_idxs <- natural_param_idxs[(length(natural_param_idxs)/2 + 1):length(natural_param_idxs)]

optim_pars[,first_natural_history_idxs] <- optim_pars[,second_natural_history_idxs]
saveRDS(optim_pars, here("inst/optims/7-19-19/top25_best_simultaneous_params.rds"))

# optim_pars[,natural_param_idxs[1:length(natural_param_idxs)/2]] <- 
# 	optim_pars[,natural_param_idxs[(length(natural_param_idxs)/2 + 1):length(natural_param_idxs)]] 

optim_pars_LA <- optim_pars[, 1:(ncol(optim_pars)/2)]
optim_pars_MA <- optim_pars[, (ncol(optim_pars)/2) + 1:(ncol(optim_pars)/2)]


colnames(optim_pars_LA) <- names(theta_la)
colnames(optim_pars_MA) <- names(theta_ma)

optim_pars_LA[, natural_history_parameters] <- optim_pars_MA[, natural_history_parameters]

# library(ggplot2)

# natural_history_parameters <- optim_pars_LA[, natural_history_parameters]

# natural_history_parameters <- reshape2::melt(natural_history_parameters)

# ggplot(natural_history_parameters, aes(x = variable, y=value, color = variable, fill = variable)) + 
#   geom_violin() + 
# 	geom_jitter() + 
# 	facet_grid(~variable) 

top5_maxs <- sort(maxs, decreasing=T)[1:5]
top5_idxs <- which(maxs %in% top5_maxs)

optim_top5_pars_LA <- optim_pars_LA[top5_idxs,]
optim_top5_pars_MA <- optim_pars_MA[top5_idxs,]

saveRDS(optim_pars_LA, file.path(files_path, 'la_25_best_pars.rds'))
saveRDS(optim_pars_MA, file.path(files_path, 'ma_25_best_pars.rds'))
saveRDS(optim_top5_pars_LA, file.path(files_path, 'la_top5_best_pars.rds'))
saveRDS(optim_top5_pars_MA, file.path(files_path, 'ma_top5_best_pars.rds'))

