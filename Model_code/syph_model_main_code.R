#####################################################
### metapopulation model of syphilis transmission ###
###            includes repeat infections         ###
###                modified for RCPP              ###
#####################################################
#setwd("~/Dropbox/PPML syphilis/")
source("Model_code/load.start.conditions.R")
state <<- "LA"   # MA=Massachusetts; LA=Louisiana
load.start()
bezier <- bezier::bezier
#theta[49] <- 0.9
#library(tmvtnorm)
#theta[28] <- log(3)
#theta[23] <- logit(0.65)
#theta[26] <- log(100)
#theta[25] <- log(50)
#test<-my_mcmc(target=dLogPosterior, init.theta=theta,proposal.sd=sd.theta/100,n.iterations=10000,limits=list(low, hi))
#test<-my_mcmc(target=dLogPosterior, init.theta=theta,proposal.sd=sd.theta/100,n.iterations=10000)
#save(test, file=paste("Model_outputs/test_mcmc_MA_constrain_dur_imm_early_2_5", Sys.time()))
#browser()
###################################################################################################
### to load outputs from MCMC and visualize calibration targets and parameter priors/posteriors ###
###################################################################################################

###change this to the apppropriate file name (stored in 'Model_outputs' folder)
#file.out <- "syph_test_p_early_only_LA_2018-03-07 21-48-35"
file.out <- "test_mcmc_LA_constrain_dur_imm_early_2_5_loosen_screening_msm_method6 2018-09-26 20-36-19"
#file.out <- "test_mcmc_MA_constrain_dur_imm_early_2_5 2018-12-14 20-24-46"

#check to make sure data are loaded from the correct site for data visualization
site.new <- ifelse(grepl("MA", file.out), "MA", "LA")
if(site.new!=state) {
  state <- site.new
  load.start()
  bezier <- bezier::bezier 
}

#load the trace, check the results, and burn and trim the MCMC chain
source("Model_code/mcmc.fit.function.R") #sample from posterior distributions to compare model outputs to calibration targets 
source("Model_code/syph_model_outputs.R") # visualize model outputs after running calibration and save to pdf
load(paste("Model_outputs/", file.out, sep="") ) #read in saved mcmc outputs
my_trace<-coda::mcmc(test$trace)
1-coda::rejectionRate(my_trace) #acceptance rate
xyplot(my_trace[,1:10])

#plotESSBurn(my_trace)
#burn.val <- max(raftery.diag(my_trace)$resmatrix[,1]) #use raftery diagnostic to get required burn-in for traces, pick the max for safety
burn.val <- 500 #2000
trace.burn <- fitR::burnAndThin(my_trace,burn=burn.val)
trace.burn.thin <- burnAndThin(trace.burn, thin=2) #original value 60
trace<-trace.burn.thin

post.sample<-model.fits(trace,100)
#browser()
plot.posteriors(post.sample)
#theta <- unlist(post.sample$theta.list[5,])
#save(theta, file="R_inputs/theta_LA.rda")

###################
### PREDICTIONS ###
###################
############################################################
### run intervention scenarios with posterior parameters ###
############################################################

#run the model comparing different interventions
#Current interventions
#1: screening at 2016 levels
#2: screen MSM at guideline levels, rest of population at 2016 levels
#3: screen MSM at guideline levels, rest of population annually
#4: Screen entire population annually
#5: Screen prior infecteds every 3 months, rest of population at 2016 levels
#6: Screen prior infecteds every 3 months, rest of population annually


#setwd("~/Dropbox/PPML syphilis/")
#source("Model_code/load.start.conditions.interv.R")
#state <- "LA"   # MA=Massachusetts; LA=Louisiana
#load.start()
#bezier <- bezier::bezier 

#n.interv = 6  #number of interventions

# ### uncomment to run model for the 6 intervention scenarios with samples drawn from parameter posteriors ###
#load("Model_outputs/Main/merged_trace")
#post.sample<-model.runs(trace, sample.size=2) 
#save(post.sample, file="Model_outputs/Main/pred_interv")

### load saved model runs to plot ###
#load("Model_outputs/Main/pred_interv")  #load in saved outputs for different screening interventions
#plot.interv()
#plot.interv.comp()  # this has not been generated yet !!! (compare intervention effectiveness?)






### plot immunity parameters ###
#png(file=paste("Model_outputs/dur_imm_",state, ".png", sep=""), units="in", width=10, height=5, res=300)
#grid.arrange(plot.dur.imm.inf, plot.dur.imm.early, plot.dur.immune,
             #ncol=3)
#dev.off()
##############################################