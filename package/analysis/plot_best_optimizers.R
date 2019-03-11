library(syphLAMA)

state <<- 'LA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', 'la_top5_pars.rds', package='syphLAMA'))
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir='~/Desktop/')


state <<- 'MA'
load.start()
trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', 'ma_top5_pars.rds', package='syphLAMA'))
post.sample <- model.fits(trace, use_trace_without_sampling=T)
plot.posteriors(post.sample, output_dir='~/Desktop/')
