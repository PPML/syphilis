# Run AMCMC With Top Values from Optimization
library(adaptMCMC)
library(here)
devtools::load_all()

args <- commandArgs(trailingOnly = T)
array_id <- args[[1]]
output_dir <- here("inst/mcmc/9-24-19/")


### Construct dLogPosterior_simultaneous to calibrate to both 
source(here("analysis/optimize/construct_simultaneous_optimization_function.R"))


top_optim_trace_path_la <- paste0(here(paste0("inst/optims/9-23-19/la_top5_best_pars.rds")))
optim_trace_la <- readRDS(top_optim_trace_path_la)
theta_la <- unlist(optim_trace_la[ (as.numeric(array_id) %% 5)+1, ])

top_optim_trace_path_ma <- paste0(here(paste0("inst/optims/9-23-19/ma_top5_best_pars.rds")))
optim_trace_ma <- readRDS(top_optim_trace_path_ma)
theta_ma <- unlist(optim_trace_ma[ (as.numeric(array_id) %% 5)+1, ])

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
