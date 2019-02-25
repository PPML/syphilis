###########################################################
###     function to sample from mcmc parameter sets     ###
###########################################################

#to compare model outputs to calibration targets
#trace = burned and thinned mcmc trace
#sample size = number of samples to draw (want 1000+ for actual analysis)

model.fits<- function(trace, sample.size){
  sample.size <- min(c(sample.size, nrow(trace)))
  ind <- sample(1:nrow(trace), sample.size, replace = TRUE)
  names(ind) <- ind
  outputs = NULL
  theta.list = NULL
  for(i.iter in 1:length(ind)) {
    theta.list <- rbind(theta.list, trace[ind[i.iter],colnames(trace)])
    #out_all <- rbind(out_all, outputs[ind[i.iter],])
  }
  for(i.iter in 1:nrow(theta.list)) {
    res<-model.epi.loglik(unlist(theta.list[i.iter,]))
    outputs <- rbind(outputs, c(i.iter,res))
  }
  #browser()
  return(list(theta.list=theta.list, outputs=outputs))
}