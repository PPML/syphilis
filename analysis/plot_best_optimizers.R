devtools::load_all()
library(here)

state <<- 'LA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', '6-6-19', 'la_25_best_pars.rds', package='syphilis'))
post.sample <- model.fits(trace, use_trace_without_sampling=T)
showCounterfactual <- FALSE
plot.posteriors(post.sample, output_dir=here('inst/optims/'))


state <<- 'MA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', '6-6-19', 'ma_25_best_pars.rds', package='syphilis'))
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir=here('inst/optims/'))
