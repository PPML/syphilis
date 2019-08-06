devtools::load_all()

# Since over-weighting the %Diagnoses among MSM and %Diagnoses which are Early
# likelihoods doesn't seem to be working, I'm going to try optimizing just the
# case report trend likelihoods and the MSM likelihood here to see how that
# turns out without the constraints of satisfying the other likelihoods.

test_likelihood <- function(theta) { 
  pred <- prediction.epi(theta)

  beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM

  ll.msm <- sum(dbeta(x=p.msm.cases,beta.params.msm$alpha/10, beta.params.msm$beta/10, log=TRUE)) #log likelihood for proportion of male cases in MSM

  print(ll.msm)
  return(ll.msm)
}

test_dLogPosterior <- function(theta) { 

  if (length(names(theta)) == 0) stop("theta must have names.")
  tryCatch({
		log.prior <- sum(prior_components(theta))
		log.likelihood <- test_likelihood(theta)
		log.posterior <- log.prior + log.likelihood
    print(log.posterior)
		if (is.na(log.posterior)||is.nan(log.posterior)) log.posterior <- -1e32
			return(as.numeric(log.posterior))
		},
  error = function(x) -1e32)
}


state <- 'MA'
load.start()

dLogPosterior(theta_ma)
test_likelihood(theta_la)
test_likelihood(theta_ma)


  out <- optim(
    par = theta_ma,
    fn = dLogPosterior,
    method = "Nelder-Mead",
    control = list(fnscale = -1, maxit=1000)
  )


new_theta <- out$par
names(new_theta) <- names(theta_ma)

saveRDS(new_theta, "new_theta.rds")

trace <- rbind.data.frame(new_theta, new_theta, new_theta)
colnames(trace) <- names(theta_ma)

trace.burn.thin <- trace.burn <- trace
post.sample <- model.fits(trace, use_trace_without_sampling=T)

plot.posteriors(post.sample, "./")
