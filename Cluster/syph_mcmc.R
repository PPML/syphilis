### syphilis model ###
setwd("~/Rfiles/syph")
source("Functions/load.start.conditions.rep.R")
state <- "LA"
load.start()
bezier <- bezier::bezier

test<-my_mcmc(dLogPosterior, theta,sd.theta/100, 10000)
save(test, file=paste("syph_test_all_",state,"_", Sys.time(), sep=""))
