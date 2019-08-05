devtools::load_all()
library(here)

state <<- 'LA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', '7-19-19', 'la_25_best_pars.rds', package='syphilis'))
trace <- rbind.data.frame(theta_la, theta_la, theta_la, theta_la, theta_la)
colnames(trace) <- names(theta_la)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
showCounterfactual <- FALSE
plot.posteriors(post.sample, output_dir="~/Desktop/")


state <<- 'MA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', '7-19-19', 'ma_25_best_pars.rds', package='syphilis'))
trace <- rbind.data.frame(theta_ma, theta_ma, theta_ma, theta_ma, theta_ma)
colnames(trace) <- names(theta_ma)
trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir="~/Desktop/")
