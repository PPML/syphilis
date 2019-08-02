
          ############################
          #          Setup           #
          ############################

# Dependencies

devtools::load_all()
library(here)

# Capture command line arguments 
args <- commandArgs(trailingOnly=TRUE)

# Get the array_id from SBATCH as nth_sim
nth_sim <- as.numeric(commandArgs(trailingOnly=T)[[1]])
if (is.null(nth_sim)) nth_sim <- 1

# For example, the script might be run as 
# Rscript optimize_with_overweighted_likelihoods.R overweight_msm_lik

# OR
# Rscript optimize_with_overweighted_likelihoods.R overweight_early_lik

# OR
# Rscript optimize_with_overweighted_likelihoods.R overweight_msm_lik overweight_early_lik

option1 <- args[[2]]
option1 <- as.character(option1)
assign(option1, TRUE)

if (length(args) >= 3) {
  option2 <- args[[3]]
  assign(option2, TRUE)
}

# script config
N_loops <- 20
output_dir_base <- "analysis/optimize_experiment/optims/"
output_dir_extension <- 
  if (exists("overweight_msm_lik") & exists("overweight_early_lik")) { 
    "both/"
  } else if (exists("overweight_msm_lik")) { 
    "overweight_msm_lik/"
  } else if (exists("overweight_early_lik")) { 
    "overweight_early_lik/"
  } else {
    stop("Neither overweight_msm_lik nor overweight_early_lik exist!")
  }

output_directory <- here(paste0(output_dir_base, output_dir_extension))

         #############################################
         #    construct dLogPosterior_simultaneous   # 
         #############################################

# Construct the simultaneous optimization function
source(here("analysis/optimize/construct_simultaneous_optimization_function.R"))


             #################################
             #     ~  o p t i m i z e  ~     #
             #################################

# main optim loop -- run optims until 'done' and checkpoint outputs iteratively 
n_done <- 0
done <- F
optim_list <- list()
current_val <- dLogPosterior_simultaneous(theta_simultaneous)

while (!done) {
  out <- optim(
    par = theta_simultaneous,
    fn = dLogPosterior_simultaneous,
    method = "Nelder-Mead",
    control = list(fnscale = -1, maxit=1000)
  )

  if (out$value > -1e32 && out$value > current_val) {
		# update and save last optim loop
    n_done <- n_done + 1
    done = n_done == N_loops
    optim_list[[n_done]] <- out
    theta <- out$par
		current_val <- out$value
    saveRDS(optim_list,
            paste0(
              output_directory,
              "optim_list_",
              nth_sim,
              ".rds"
            )
    )
  }
}


