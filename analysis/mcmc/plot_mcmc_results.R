library(here)
devtools::load_all()

output_directory <- here("inst/mcmc/5-2-19/")
for (state in c("LA", "MA")) {
  load.start()
  trace <- trace.burn.thin <- trace.burn <- readRDS(
    file.path(output_directory, paste0("trace_burn_and_thinned_", state, ".rds")))
  
  post.sample <- model.fits(trace, use_trace_without_sampling = T)
  pred <- as.data.frame(post.sample$outputs)
  saveRDS(post.sample, paste0(output_directory, state, "_posterior_sample.rds"))
  plot.posteriors(output_directory, post.sample = post.sample)
}
