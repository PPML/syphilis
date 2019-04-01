library(syphLAMA)

# Find optimization.R output files
files_path <- '~/2019/February/26/optim/'
la_files <- grep('LA', list.files(files_path, full.names=T), value=T)
ma_files <- grep('MA', list.files(files_path, full.names=T), value=T)

# Read them
la_optims <- lapply(la_files, readRDS)
ma_optims <- lapply(ma_files, readRDS)

# Write helpers to extract data
get_optim_vals <- function(optim_chain) {
  sapply(optim_chain, `[[`, 'value')
}

get_optim_pars <- function(optim_chain) {
  lapply(optim_chain, `[[`, 'par')
}

# Extract posterior values
la_optim_vals <- lapply(la_optims, get_optim_vals)
ma_optim_vals <- lapply(ma_optims, get_optim_vals)


# Find the best posterior values
la_maxs <- sapply(la_optim_vals, max) 
la_max_idxs <- sapply(la_optim_vals, function(x) which(x==max(x))[[1]])

ma_maxs <- sapply(ma_optim_vals, max) 
ma_max_idxs <- sapply(ma_optim_vals, function(x) which(x==max(x))[[1]])

print(sort(la_maxs, decreasing=T))
print(sort(ma_maxs, decreasing=T))

# Get the corresponding parameter vectors for the top 5 posterior values
# get_optim_pars -> which are in top 5 -> 
la_top5_maxs <- sort(la_maxs, decreasing=T)[1:5]
ma_top5_maxs <- sort(ma_maxs, decreasing=T)[1:5]

# which chains are best 5 
la_top5_idxs <- which(la_maxs %in% la_top5_maxs)
ma_top5_idxs <- which(ma_maxs %in% ma_top5_maxs)

# get each chain's param list
la_optim_pars <- lapply(la_optims, get_optim_pars)
ma_optim_pars <- lapply(ma_optims, get_optim_pars)

# get each chain's best param
la_optim_max_pars <- lapply(1:length(la_optim_pars), function(x) la_optim_pars[[x]][[la_max_idxs[[x]]]])
ma_optim_max_pars <- lapply(1:length(ma_optim_pars), function(x) ma_optim_pars[[x]][[ma_max_idxs[[x]]]])

# format every chains best params as dataframe
param_names <- names(la_optim_max_pars[[1]]) 
la_optim_max_pars <- unname(do.call(rbind.data.frame, la_optim_max_pars)) 
colnames(la_optim_max_pars) <- param_names
ma_optim_max_pars <- unname(do.call(rbind.data.frame, ma_optim_max_pars)) 
colnames(ma_optim_max_pars) <- param_names

# get top 5
la_top5_pars <- la_optim_max_pars[la_top5_idxs,]
ma_top5_pars <- ma_optim_max_pars[ma_top5_idxs,]


# Save results
outdir <- "~/2019/March/11/"
dir.create(outdir, recursive=T)
saveRDS(la_top5_pars, file.path(outdir, 'la_top5_pars.rds'))
saveRDS(ma_top5_pars, file.path(outdir, 'ma_top5_pars.rds'))

