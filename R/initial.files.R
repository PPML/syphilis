####################################################
### read in required data files and prepare data ###
####################################################

# all data files are stored in 'R_inputs folder ###
# refer to 'file_descriptions.docx' in 'R_inputs' folder for details about files
# state can be either "MA" (Massachusetts) or "LA" (Louisiana)

#' Load Initial Files by State (LA/MA)
initial.files <- function(state) {
    age.dist.dat <<- read.delim(system.file("extdata/age_nsfg.txt", package='syphilis')) #proportion of partnerships with person of same age group
    diag.age.sex.dat <<-as.matrix(read.delim(system.file(paste("extdata/diag_age_sex_rate_", state, ".txt", sep=""), package='syphilis'), row.names=1)) #case rates by age sex only, 2000-2015, per 100,000 population
    diag.rr <<- read.delim(system.file(paste("extdata/case_subpop_rr_", state, ".txt", sep=""),package='syphilis')) ##RR of diagnosis rate for black and Hispanic, relative to overall population of same age and sex (5-yr avg, 2012-2016)
    msm.dat <<- as.matrix(read.delim(system.file(paste("extdata/msm_hiv_proportions_", state, ".txt", sep=""), package='syphilis'), row.names=1)) #proportion of males cases in MSM, and proportion of HIV+ among MSM cases, 2007-2016
    stage.dat <<- as.matrix(read.delim(system.file(paste("extdata/stage_diag_", state, ".txt", sep=""), package='syphilis'), row.names=1)) #proportion of cases diagnosed as secondary or early latent, 2000-2016, excluding F 45+ due to small case numbers
    p.early.dat <<- read.delim(system.file(paste("extdata/p_early_", state, ".txt", sep=""), package='syphilis')) #proportion of diagnosed syphilis that are early (of all cases, including late latent, 2011-2015 for all cases in given state from cdc)
    diag.subpop <<- tibble::rownames_to_column(as.data.frame(read.delim(system.file(paste("extdata/subpop_diag_", state, ".txt", sep=""), package='syphilis'), row.names=1)), "year" )#diagnosis rates by subpop

    load(system.file(paste("extdata/theta_", state, ".rda", sep=""),package='syphilis')) #load starting values of parameters
    theta <<- theta
    
    #prep data used for calibration
    diag.age.sex.rate.m <<- as.vector(diag.age.sex.dat[, c("diag_all_y_m", "diag_all_o_m")])  #diagnosed early syphilis  rates per 100,000, for males 
    diag.age.sex.rate.m.sd <<- as.vector(diag.age.sex.dat[, c("diag_all_y_m_sd", "diag_all_o_m_sd")])  #SD for diagnosed case data for males 
    diag.age.sex.rate.m.sd <<- diag.age.sex.rate.m.sd * 0.25
    diag.age.sex.rate.f.y <<- as.vector(diag.age.sex.dat[, "diag_all_y_f"])  #diagnosed early syphilis  rates per 100,000, for young females 
    diag.age.sex.rate.f.y.sd <<- as.vector(diag.age.sex.dat[, "diag_all_y_f_sd"])  #SD for diagnosed case data for young females 
    diag.age.sex.rate.f.o <<- as.vector(diag.age.sex.dat[, "diag_all_o_f"])  #diagnosed early syphilis  rates per 100,000, for older females 
    diag.age.sex.rate.f.o.sd <<- as.vector(diag.age.sex.dat[, "diag_all_o_f_sd"])  #SD for diagnosed case data for older females 
    rr.diag.subpop <<- as.numeric(diag.rr[,"mean"]) 
    rr.diag.subpop.sd <<- as.numeric(diag.rr[,"sd"])
    p.msm.cases <<- as.vector(msm.dat[,c("pMSM_y","pMSM_o")]) #proportion of MALE cases diagnosed in MSM
    p.msm.cases.var <<- as.vector(msm.dat[,c("pMSM_y_var", "pMSM_o_var")]) #variance 
    p.hiv.cases <<- as.vector(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"])))) #proportion of MALE cases diagnosed in MSM
    p.hiv.cases.var <<- as.vector(subset(msm.dat[,c("pHIV_y_var","pHIV_o_var")], (!is.na(msm.dat[,"pHIV_y_var"])) & (!is.na(msm.dat[,"pHIV_o_var"])))) #variance 
    p.sec <<- as.vector(stage.dat[,c("pSec_y_m","pSec_o_m","pSec_y_f")]) #proportion cases diagnosed with secondary syphilis
    p.sec.var <<- as.vector(stage.dat[,c("pSec_var_y_m","pSec_var_o_m","pSec_var_y_f")]) 
    p.el <<- as.vector(stage.dat[,c("pEL_y_m","pEL_o_m","pEL_y_f")]) #proportion cases diagnosed with early latent syphilis
    p.el.var <<- as.vector(stage.dat[,c("pEL_var_y_m","pEL_var_o_m","pEL_var_y_f")]) 
    #browser() 
    #read in model priors (used by model.priors.fun.R)
    priors <<- read.delim(system.file("extdata/priors.txt", package='syphilis'))  # read in priors 
    sd.theta <<- priors$sd.transf.1 #starting value for standard deviation associated with each parameter, adapted during fitting
    prior.param1 <<-  priors$param1 #first parameter describing probablity distributions
    names(prior.param1) <<-priors$parameter
    prior.param2 <<-  priors$param2 #second parameter describing probability distributions
    names(prior.param2) <<- priors$parameter
  }
