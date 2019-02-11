#####################################################
### metapopulation model of syphilis transmission ###
###            includes repeat infections         ###
###                modified for RCPP              ###
#####################################################
#setwd("~/Dropbox/PPML syphilis/")
source("Model_code/load.start.conditions.interv.R")
state <- "LA"   # MA=Massachusetts; LA=Louisiana
load.start()
bezier <- bezier::bezier 

#test<-my_mcmc(target=dLogPosterior, init.theta=theta,proposal.sd=sd.theta/100,n.iterations=100)
#save(test, file=paste("test_mcmc", Sys.time()))

###################################################################################################
### to load outputs from MCMC and visualize calibration targets and parameter priors/posteriors ###
###################################################################################################

###change this to the apppropriate file name (stored in 'Model_outputs' folder)
file.out <- "test_mcmc_LA_constrain_dur_imm_early_2_5_loosen_screening_msm_f_method6 2018-09-28 15-13-26"

#check to make sure data are loaded from the correct site for data visualization
site.new <- ifelse(grepl("MA", file.out), "MA", "LA")
if(site.new!=state) {
  state <- site.new
  load.start()
  bezier <- bezier::bezier 
}
n.interv = 6  #number of interventions
#load the trace, check the results, and burn and trim the MCMC chain
#source("Model_code/mcmc.fit.function.R") #sample from posterior distributions to compare model outputs to calibration targets 
#source("Model_code/syph_model_outputs.R") # visualize model outputs after running calibration and save to pdf
source("Model_code/mcmc.fit.interv.fun.R") #function to iterate through different interventions and print outputs
source("Model_code/plot.interv.R") #plot comparison of screening interventions
source("Model_code/model.out.interv.fun.R")
source("Model_code/model.pred.interv.fun.R")

load(paste("Model_outputs/", file.out, sep="") ) #read in saved mcmc outputs


my_trace<-coda::mcmc(test$trace)
#trace2<- my_trace
#my_trace<-mcmc(test$trace)
1-coda::rejectionRate(my_trace) #acceptance rate
xyplot(my_trace[,1:10])

#plotESSBurn(my_trace)
trace.burn <- burnAndThin(my_trace,burn=2000)
#xyplot(x=trace.burn[,1:15])
#effectiveSize(trace.burn)
#acfplot(x=trace.burn, lag.max=60)
trace.burn.thin <- burnAndThin(trace.burn, thin=60)
trace<-trace.burn.thin
#xyplot(x=trace.burn.thin[,1:15])
#effectiveSize(trace.burn.thin)
#acfplot(x=trace.burn.thin, lag.max=60)
#summary(trace.burn)
#summary(trace.burn.thin)
#plotPosteriorDensity(list(unthinned=trace.burn, thinned=trace.burn.thin))

#out <-test$outputs #model outputs from MCMC


#post.sample<-model.runs(trace,100)
post.sample<-model.runs(trace,10)
plot.interv.comp()
#theta <- unlist(post.sample$theta.list[5,])
#save(theta, file="R_inputs/theta_LA.rda")

######

pdf(file=paste("MCMC_outputs/dur_imm_",state, ".pdf", sep=""), width=8, height=5, paper="USr")

png(file=paste("MCMC_outputs/dur_imm_",state, ".png", sep=""), units="in", width=10, height=5, res=300)
grid.arrange(plot.dur.imm.inf, plot.dur.imm.early, plot.dur.immune,
             ncol=3)
dev.off()
#save(my_trace, file="mcmc_28_Jun_2016")
load("MCMC_outputs/mcmc_prev_pmsm_diagearlyall01pct_lateFonly10pct_20K 2016-11-14 03:26:31")
#to combine traces:
trace3<-mcmc(test$trace)
# trace1<- mcmc(traceObj1$trace)
# trace2<- mcmc(traceObj2$trace)
trace <- mcmc.list(list(trace1,trace2,trace3))
1-rejectionRate(trace)
effectiveSize(trace)
trace.burn <- burnAndThin(trace, burn = 5000)
effectiveSize(trace.burn)
acfplot(trace.burn, lag.max=60)
trace.burn.thin<-burnAndThin(trace.burn, thin=60)
effectiveSize(trace.burn.thin)
#xyplot(trace.burn.thin)
plotPosteriorDensity(list(chain1 = trace.burn.thin[[1]][,1:25], chain2 = trace.burn.thin[[2]][,1:25], chain3=trace.burn.thin[[3]][,1:25]))
plotPosteriorDensity(list(chain1 = trace.burn.thin[[1]][,26:50], chain2 = trace.burn.thin[[2]][,26:50], chain3=trace.burn.thin[[3]][,26:50]))
plotPosteriorDensity(list(chain1 = trace.burn.thin[[1]][,51:75], chain2 = trace.burn.thin[[2]][,51:75], chain3=trace.burn.thin[[3]][,51:75]))
plotPosteriorDensity(list(chain1 = trace.burn.thin[[1]][,76:99], chain2 = trace.burn.thin[[2]][,76:99], chain3=trace.burn.thin[[3]][,76:99]))
gelman.diag(trace.burn.thin) #want point est <1.1
gelman.plot(trace.burn.thin)

densityplot(trace.burn.thin)
levelplot(trace.burn.thin[[1]],col.regions=heat.colors(100))

### check correlations
levelplot(trace.burn.thin[,62:75],col.regions=heat.colors(100))
cor(trace.burn.thin[,62:75])
#PRCC
require(sensitivity)
pred.theta<-as.data.frame(post.sample$theta)
y<-pred[,"prev.fit.prev.extra13"]
x<-pcc(pred.theta, y, rank=TRUE)
