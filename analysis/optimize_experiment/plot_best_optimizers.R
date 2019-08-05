devtools::load_all()
library(here)

files_path <- here("analysis/optimize_experiment/optims/overweight_msm_lik/")

state <<- 'LA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(paste0(files_path, 'la_25_best_pars.rds'))
trace <- rbind.data.frame(theta_la, theta_la, theta_la, theta_la, theta_la)
colnames(trace) <- names(theta_la)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
showCounterfactual <- FALSE
plot.posteriors(post.sample, output_dir=files_path)


state <<- 'MA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(paste0(files_path, 'ma_25_best_pars.rds'))
trace <- rbind.data.frame(theta_ma, theta_ma, theta_ma, theta_ma, theta_ma)
colnames(trace) <- names(theta_ma)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir=files_path)
