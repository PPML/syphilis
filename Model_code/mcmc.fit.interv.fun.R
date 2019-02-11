###########################################################
###     function to sample from mcmc parameter sets    ###
###########################################################

#to compare model outputs under different screening scenarios
#trace = burned and thinned mcmc trace
#sample size = number of samples to draw (want 1000+ for actual analysis)
#loops through the different interventions
  
  model.runs<-function(trace, sample.size){
    sample.size <- min(c(sample.size, nrow(trace)))
    ind <- sample(1:nrow(trace), sample.size, replace = TRUE)
    names(ind) <- ind
    outputs = NULL
    theta.list = NULL
    for(i.iter in 1:length(ind)) {
      theta.list <- rbind(theta.list, trace[ind[i.iter],colnames(trace)])
    }
    for(i.iter in 1:nrow(theta.list)) {
      for(i.interv in 1:n.interv) {
        interv=i.interv
        res<-prediction.epi(unlist(theta.list[i.iter,]), interv)
        outputs <- rbind(outputs, c(iter=i.iter,res, interv=interv))
      }
    }
    
    return(list(theta.list=theta.list, outputs=outputs))
  }
  