devtools::load_all()
library(here)

optim_path <- here("inst/optims/10-12-19/")

state <<- 'LA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(paste0(optim_path, 'la_top5_best_pars.rds'))
colnames(trace) <- names(theta_la)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
showCounterfactual <- FALSE
plot.posteriors(post.sample, output_dir=optim_path)


state <<- 'MA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(paste0(optim_path, 'ma_top5_best_pars.rds'))
colnames(trace) <- names(theta_ma)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir=optim_path)
