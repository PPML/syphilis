library(here)
library(coda)
library(fitR)

output_directory <- here("inst/mcmc/5-2-19/")

for (state in c("LA", "MA")) {
  calibration_files <- grep(paste0(state, "_mcmc"), list.files(output_directory, full.names = T), value = T)
  mcmc_list <- mcmc.list(lapply(calibration_files, function(f) mcmc(readRDS(f)$samples)))
  mcmc_list <- burnAndThin(mcmc_list, burn = 15000, thin = 100)
  trace <- do.call(rbind.data.frame, mcmc_list)
  trace <- dplyr::sample_n(trace, 1000)
  trace.burn.thin <- trace.burn <- trace
  saveRDS(trace.burn.thin, file.path(output_directory, paste0("trace_burn_and_thinned_", state, ".rds")))
}
