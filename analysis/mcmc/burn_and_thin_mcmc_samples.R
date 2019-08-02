library(here)
library(coda)
library(fitR)
devtools::load_all()

state <- 'LA' 
# Choice of LA or MA doesn't matter here, we only load.start() so we have 
# access to prior.param1 and prior.param2 which don't vary by state
load.start()


output_directory <- here("inst/mcmc/8-1-19/")

calibration_files <- grep(paste0("mcmc_"), list.files(output_directory, full.names = T), value = T)
mcmc_list <- mcmc.list(lapply(calibration_files, function(f) mcmc(readRDS(f)$samples)))
# mcmc_list <- burnAndThin(mcmc_list, burn = 15000, thin = 100)
trace <- do.call(rbind.data.frame, mcmc_list)


trace_la <- trace[,1:(ncol(trace)/2)]
trace_ma <- trace[,(ncol(trace)/2) + 1:(ncol(trace)/2)]

trace_la[,natural_history_parameters] <- trace_ma[,natural_history_parameters]


within <- function(x, y, z) { return(x <= z & x >= y) } 

valid_rows <- 
within(exp(trace_la["log.dur.incub"]),prior.param1["dur.incub"],prior.param2["dur.incub"]) &
within(exp(trace_la["log.dur.prim"]),prior.param1["dur.prim"],prior.param2["dur.prim"]) &
within(exp(trace_la["log.dur.sec"]),prior.param1["dur.sec"],prior.param2["dur.sec"]) &
within(exp(trace_la["log.dur.imm.early"]),prior.param1["dur.imm.early"],prior.param2["dur.imm.early"])

trace_la <- trace_la[valid_rows,]
trace_ma <- trace_ma[valid_rows,]


sample_idxs <- sample.int(nrow(trace_la), 1000)
trace_la <- trace_la[sample_idxs, ]
trace_ma <- trace_ma[sample_idxs, ]


# saveRDS(trace.burn.thin, file.path(output_directory, "simultaneous_trace_burn_and_thinned.rds"))

usethis::use_data(trace_la, overwrite=TRUE)
usethis::use_data(trace_ma, overwrite=TRUE)

saveRDS(trace_la, here("inst/mcmc/8-1-19/trace_la_burned_and_thinned.rds"))
saveRDS(trace_ma, here("inst/mcmc/8-1-19/trace_ma_burned_and_thinned.rds"))
