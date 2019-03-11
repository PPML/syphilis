library(syphLAMA)

# command line arguments
site <<- commandArgs(trailingOnly=T)[[1]]
state <<- site
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[2]])

# script config
N_loops <- 20
output_directory <- "~/2019/February/26/optim/"

# setup syphLAMA load trace
load.start()
orig_theta_names <- names(theta)
x <- load(system.file(paste0("mcmc/mcmc_", site), package='syphLAMA'))
post.sample <- get(x = x); rm(x)
trace <- post.sample$trace
theta <- unlist(trace[sample.int(nrow(trace), 1),])

# make sure the theta vector is named properly 
dLogPosterior_for_optim <- function(theta) {
	names(theta) <- orig_theta_names
	d <- dLogPosterior(theta)
	print(d)
	return(d)
}

# main optim loop -- run optims until 'done' and checkpoint outputs iteratively 
n_done <- 0
done <- F
optim_list <- list()
while (!done) {
  out <- optim(
    par = theta,
    fn = dLogPosterior_for_optim,
    method = c("BFGS", "Nelder-Mead")[[if (n_done > 5) 1 else 2]],
    control = list(fnscale = -1, maxit=1000)
  )

  if (out$value != -1e32) {
		# update and save last optim loop
    n_done <- n_done + 1
    done = n_done == N_loops
    optim_list[[n_done]] <- out
    theta <- out$par
    saveRDS(optim_list,
            paste0(
              output_directory,
              site,
              "_optim_list_",
              nth_sim,
              ".rds"
            )
    )
  } else {
    theta <-  unlist(trace[sample.int(nrow(trace), 1),])
 }
}
