# Adapted from analysis/plot_best_optimizers.R

# I want to determine whether or not the duration of immunity is causing 
# some of the jaggedness that we've seen in the optimized model fits.

# To do this, we'll first generate some plots demonstrating the jagged 
# trends we've observed in the outputs plot generated from the optimized
# parameter vectors. 
# 
# Then we will set the durations of immunity to very low values and 
# observe how the simulated trajectories change.

devtools::load_all()

state <<- 'LA'
load.start()
trace <- readRDS(system.file('optims', 'la_top5_pars.rds', package='syphLAMA'))
theta <- unlist(trace[1,])
# sol <- 

# post.sample <- model.fits(trace, use_trace_without_sampling=T)

# plot.posteriors(post.sample, output_dir='~/Desktop/')

# state <<- 'MA'
# load.start()
# trace.burn.thin <- trace.burn <- trace <- readRDS(system.file('optims', 'ma_top5_pars.rds', package='syphLAMA'))
# post.sample <- model.fits(trace, use_trace_without_sampling=T)
# plot.posteriors(post.sample, output_dir='~/Desktop/')
