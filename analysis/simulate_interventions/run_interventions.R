devtools::load_all(".")
library(here)

output_path <- here('inst/interventions/')

state <- commandArgs(trailingOnly=TRUE)[[1]]

state <- toupper(state)

load.start()

trace <- get(paste0('trace_', tolower(state)))

df <- simulate_interventions_trace(trace, info_extractor = extract_outcomes)

saveRDS(df, file = file.path(output_path, paste0(tolower(state), '_interventions_df.rds')))

