##############################################
###   additional functions used for MCMC  ####
###   and other miscellaneous functions   ####
##############################################

#calculate likelihood of given iteration and update parameter set used
likelihood <- function(theta){ 
  Res <-  model.epi.loglik(theta)
  LogLL <- Res["ll"] #Model is returning the Log likelihood
  theta.last<<-theta
  return(list( LogLL=LogLL))
}

#calculate posterior likelihood
dLogPosterior <-function(theta) {
  if (length(names(theta)) == 0) stop("theta must have names.")
  tryCatch({
		log.prior <- sum(prior_components(theta))
		log.likelihood <- unlist(likelihood(theta)$LogLL)
		log.posterior <- log.prior + log.likelihood
		if (is.na(log.posterior)||is.nan(log.posterior)) log.posterior <- -1e32
			return(as.numeric(log.posterior))
		},
  error = function(x) -1e32)
}

#run adaptive Metroposlis Hastings MCMC algorithm (adapted from fitR package)
#see mcmcMH in fitR package for details 
my_mcmc<-function(target, init.theta, proposal.sd, n.iterations) {
  out_all <<- NULL
  trace <- myMCMC(target = target,
                  init.theta = init.theta,
                  proposal.sd = proposal.sd,
                  n.iterations = n.iterations,
                  adapt.size.start = 100,
                  adapt.shape.start = 200,
                  adapt.size.cooling=0.999,
                  verbose=FALSE)
  return(trace)
}

# function to estimate parameters for beta distribution, where pred.prop = probability, pop = population size, and scale is a scaling factor that determines variance
estBetaParams <- function(pred.prob, pop, scale) {
  alpha <- pred.prob*pop/scale
  beta <- (1-pred.prob)*pop/scale
  return(beta.params = list(alpha = alpha, beta = beta))
}

estBetaParamsVar <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  list(alpha = alpha, beta = beta)
}


logit<-function(x) {log(x/(1-x))}

ilogit <-function(x) {1/(1+exp(-x))}

# function to calculate annual rates from cumulative output, when generating annual outputs
annual.rates <- function(x){  
  annual <-rbind(rep(0,ncol(x)), apply(x,2,diff))
}

# untransform_theta
untransform_theta <- function(theta) { 
  for (iter in 1:length(theta)) { 
    if (grepl("logit\\.", names(theta)[[iter]])) { 
      # untransform logit 
      theta[[iter]] <- ilogit(theta[[iter]])

      # rename
      temp_name <- strsplit(names(theta)[[iter]], "\\.")[[1]]
      temp_name <- paste0(temp_name[2:length(temp_name)], collapse='.')
      names(theta)[[iter]] <- temp_name

    } else if (grepl("log\\.", names(theta)[[iter]])) { 
      # untransform log 
      theta[[iter]] <- exp(theta[[iter]])

      # rename
      temp_name <- strsplit(names(theta)[[iter]], "\\.")[[1]]
      temp_name <- paste0(temp_name[2:length(temp_name)], collapse='.')
      names(theta)[[iter]] <- temp_name

    } else stop("theta contains a term which is not a logit or log transformation.")
  }

  return(theta)
}
