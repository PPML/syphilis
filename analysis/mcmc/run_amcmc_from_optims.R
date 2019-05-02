# Run AMCMC With Top Values from Optimization
library(adaptMCMC)
library(here)
devtools::load_all()

args <- commandArgs(trailingOnly = T)
state <- site <- args[[1]]
array_id <- args[[2]]
output_dir <- here("inst/mcmc/5-2-19/")

load.start()

top_optim_trace_path <- paste0(here(paste0("inst/optims/", tolower(site), "_top5_pars.rds")))

optim_trace <- readRDS(top_optim_trace_path)

theta <- unlist(optim_trace[ (as.numeric(array_id) %% 5)+1, ])

sd_optims <- apply(optim_trace, 2, sd)

sample <- MCMC(
  p = dLogPosterior,
  n = 25000,
  init = theta,
  scale = sd_optims / 1000,
  adapt = T,
  acc.rate = .234, showProgressBar = T)

saveRDS(sample, file = paste0(output_dir, site, "_mcmc_25000_", array_id, ".rds"))
