#######################################################
### function to load start conditions and functions ###
#######################################################

load.start <- function(){
  library(deSolve)
  library(scales)
  library(lattice)
  library(latticeExtra)
  library(fitR)
  library(Rcpp)
  library(bezier)
  library(plyr)
  library(reshape)
  
  ### load functions ###
  source("Model_code/mixing.fun.R")  # function to calculate and balance contact matrix, outputs p.total.bal;p.total.bal.btwn;cm.low;cm.high
  source("Model_code/aging.fun.R")  #function to calculate aging parameters
  source("Model_code/population.size.fun.R") #function to calculate population sizes
  source("Model_code/model.pred.fun.R")  #function to update parameters, run model, and print outputs used for fitting
  source("Model_code/initial.files.R")  #read in required data files
  source("Model_code/bezier.funs.R") #functions to calculate time-varying parameters
  source("Model_code/myMCMC.R") #modified mcmcMH function from fitR package that ensures covariance matrix is positive definite
  source("Model_code/model.out.fun.R") #function to calculate desired model outputs for calibration
  source("Model_code/model.likelihood.fun.R") #calculate model likelihoods
  source("Model_code/mcmc.funs.R")  #functions used for mcmc algorithm
  source("Model_code/model.priors.fun.R") #calculate prior likelihood  
  sourceCpp("Model_code/syph_trans_model.cpp") #read in rcpp code for transmission model
  
  initial.files(state) #load required data files
  
  showCounterfactual <<- TRUE #turn on counterfactual comparison for contact tracting
  
  ### set up different subpopulations, sexes, and activity classes ###
  i<<-5 #number of subpopulations, 1= Black, 2=White, 3=Hispanic, 4=MSM-HIVneg, 5=MSM-HIVpos
  j<<-2 # number of sexes, 1=male, 2=female
  k<<-2 #activity class, 1=low, 2=high
  l<<-2 #age groups, 1= 20-44; 2=45-64
  index <<- i*j*k*l
  
  #### model time steps and calibration period ####
  tstep <<- 1/52  #weekly time step
  cal.period <<- 5 #duration of calibration period (2012-2016 currently)
  cal.start <<- 100  #time at which start calibration 
  model.end <<- cal.start + cal.period + 1
  end.year <<- 2016
  start.year <<- 2012
  ct.data.years <<- if(state=="MA") 12 else 6
  
  # delay to diagnosis in contact tracing
  d <<- 1/12
  
  # background antibiotic treatment
  p.abx.init <- 0.15 #set this to 0 if want to turn off abx treatment -- this implements a period (pre-cal) with high rates of abx use, to represent the intro of penicillin in the pop, which decreases to abx.background 
  p.abx.background <- 0.01
  abx.start <- cal.start-60
  abx.end <- cal.start-20
  p.abx <<- c(rep(0,abx.start), rep(p.abx.init,(abx.end-abx.start)),rep(p.abx.background, (model.end+1-abx.end)) )
  
  # when to start moving people into the prior infected compartments 
  rep.start <- cal.start - 5 #years before calibration start when start tracking prior treated infections
  rep.count <<- c(rep(0,rep.start), rep(1, model.end+1-rep.start))
  
  age.cat<<-c(25, 20) # age band widths, corresponding to 20-44 yo and 45-64 yo
  out_all <<- NULL
  omega <<- 0.5 #supply and demand of sexual partnerships, 0.5=both sexes copmromise equally, 1=females determine number of partnerhips, 0=males determine number of partnerships
  omega.t <<- if(state=="MA") matrix(c(0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0),nrow=5,byrow=T) else matrix(c(0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0),nrow=5,byrow=T)
    #parameter describing compromise in partner number across subpops --> assuming that smaller pop size determines number of partnerships, rows= subpop of F, col=subpop of M
  
  ### population characteristics for Massachusetts and Louisiana ###
  n.total <<- if(state=="MA") 4230601 else 2787423 # total population size in 2015 for (i) MA and (ii) LA
  p.msm <<- if(state=="MA") 0.045 else 0.025 # proportion of population MSM (x% of males, using estimates for MA and LA from Grey et al. 2016)
  p.s.1 <<-if(state=="MA") 0.073 else 0.32 #proportion of population in i=1 (black) in 2015
  p.s.3 <<- if(state=="MA") 0.11 else 0.05 #proportion of population in i=3 (hispanic) in 2015
  p.s.2<<-1-p.s.1-p.s.3
  p.low <<- 0.90 #proportion of population in low activity group
  p.low.msm <<-0.90 # proportion of MSM in low activity group
  p.hiv <<- if(state=="MA") c(0.05,0.12) else c(0.22, 0.28) # proportion of MSM who are HIV positive (young and old) -- 2015 estimates 
  p.sa.m.y <<- c(0.98, 0.94, 0.95, 0.94,1) #proportion of young M that are sexually active, by subpop (based on NSFG) - assuming mean for MSM and all HIV+ MSM are sexually active
  p.sa.m.o <<-c(rep(0.99,i-1),1) #proportion of M old that are sexually active (assumption, based on males aged 40-44 in NSFG reporting ever sexually active)
  p.sa.f.y <<- rep(0.96,i) #proportion of young  F that are sexually active (from NSFG)
  p.sa.f.o <<- rep(0.99,i) #proportion of old F that are sexually active (assumption, based on females aged 40-44 in NSFG reporting ever sexually active)
  
  aging.rate.m <<- ifelse(p.sa.m.y<1, (p.sa.m.o - p.sa.m.y) / (1-p.sa.m.y), 1)
  aging.rate.f <<- (p.sa.f.o - p.sa.f.y) / (1-p.sa.f.y)
  aging.rate.hiv <<- (p.hiv[2]-p.hiv[1]) / (1-p.hiv[1])
  birth.rate.hiv <<- p.hiv[1]/p.hiv[2]
  p.sa.m<<-array(c(p.sa.m.y, p.sa.m.o),dim=c(i,2)) #for each subpop, p.sa for M (y,o) 
  p.sa.f<<-array(c(p.sa.f.y, p.sa.f.o),dim=c(i,2)) #for each subpop, p.sa for M (y,o) 
  p.sa<<-c(as.vector(rep(aperm(p.sa.m, perm=c(2,1)), each=2)), as.vector(rep(aperm(p.sa.f, perm=c(2,1)), each=2)))
  
  # contact tracing probabilities and whether it is turned on
  p.ct.primsec <<- if(state == "MA") { # Black M, Hispanic M, NBNH M, Black F, Hispanic F, NBNH F
                      ct.data.ps <- read.delim("R_inputs/ct_ma_rates_ps.txt")
                      t(rbind(ct.data.ps$black.men, ct.data.ps$hispanic.men, ct.data.ps$other.men, 
                                           ct.data.ps$black.women, ct.data.ps$hispanic.women, ct.data.ps$other.women))
                  } else {
                      ct.data.ps <- read.delim("R_inputs/ct_la_rates_ps.txt")
                      t(rbind(ct.data.ps$black.men, ct.data.ps$hispanic.men, ct.data.ps$other.men, 
                                           ct.data.ps$black.women, ct.data.ps$hispanic.women, ct.data.ps$other.women))
                  }
  p.ct.el <<- if(state == "MA") { # Black M, Hispanic M, NBNH M, Black F, Hispanic F, NBNH F
                ct.data.el <- read.delim("R_inputs/ct_ma_rates_el.txt")
                t(rbind(ct.data.el$black.men, ct.data.el$hispanic.men, ct.data.el$other.men, 
                                     ct.data.el$black.women, ct.data.el$hispanic.women, ct.data.el$other.women))
              } else {
                ct.data.el <- read.delim("R_inputs/ct_la_rates_el.txt")
                t(rbind(ct.data.el$black.men, ct.data.el$hispanic.men, ct.data.el$other.men, 
                        ct.data.el$black.women, ct.data.el$hispanic.women, ct.data.el$other.women))
              }
  
  # for each subpopulation (i), need to define pop. size (M+F) and distribution of activity classes and ages
  pop.calc(n.total)
  aging.fun(age.cat)

  ### initial model conditions ###
  ind<<- as.matrix(expand.grid(1:k,1:l,1:j,1:i)) #indexing matrix (k,l,j,i)
  init.Y <<-c(rep(c(0,3,0,3),4),rep(c(0,0,0,0),1),rep(c(0,3,0,3),3),rep(c(0,0,0,0),2)) #vector of initial number of infecteds

  #initial distribution of model pop by sex, subpop, AC (M pop1 age1 L/H; M pop1 age 2 L/H; M pop2 age 1 L/H, M pop 2 age 2 L/H, etc.)
  yinit <<- c(  
    S=n.sa-init.Y,          #susceptible
    E=c(rep(0,index)),      #incubating
    I.1 = init.Y/3,         #primary syphilis
    I.2= init.Y/3,          #secondary syphilis
    L.1=init.Y/3,           #early latent syphilis
    L.2=c(rep(0,index)),    #late latent syphilis
    T1=c(rep(0,index)),     #treated, primary and secondary
    T2=c(rep(0,index)),     #treated, early latent
    T3=c(rep(0,index)),     #treated, late latent
    SR=c(rep(0,index)),     #susceptible, prior treated infection
    ER=c(rep(0,index)),     #incubating, prior treated infection
    IR.1 = c(rep(0,index)), #primary, prior treated infection
    IR.2= c(rep(0,index)),  #secondary, prior treated infection
    LR.1=c(rep(0,index)),   #early latent, prior treated infection
    LR.2=c(rep(0,index)),   #late latent, prior treated infection
    INC=c(rep(0,index)),    #cumulative incidence
    INCR = c(rep(0,index)), #cumulative incidence, prior infection
    D1=c(rep(0,index)),     #reported, primary
    D2=c(rep(0,index)),     #reported, secondary
    D3=c(rep(0,index)),     #reported, early latent
    NSA=n.nsa,              #not sexually active
    D4=c(rep(0,index)),     #reported, late latent
    DR = c(rep(0,index))    #reported, primary, secondary, early latent in prior infecteds
   )
  
  #indexes for pulling out different states/population groups from the model output matrix
  s.index <<- 1:index                      #susceptible
  e.index <<- index+1:index                #incubating
  prim.index <<- index*2+1:index           #primary infection
  sec.index <<- index*3+1:index            #secondary infection
  early.index <<- index*4+1:index          #early latent infection
  latent.index <<- index*5+1:index         #late latent infection
  treated.inf.index <<- index*6+1:index    #treated P&S
  treated.early.index <<- index*7+1:index  #treated early latent
  treated.late.index <<- index*8+1:index   #treated late latent
  sr.index <<- index*9+1:index             #re-susceptible
  er.index <<- index*10+1:index            #incubating, previously treated
  primr.index <<- index*11+1:index         #primary infection, previously treated
  secr.index <<- index*12+1:index          #secondary infection, previously treated
  earlyr.index <<- index*13+1:index        #early latent infection, previously treated
  latentr.index <<- index*14+1:index       #late latent infection, previously treated
  inc.index <<- index*15+1:index           #incidence, all
  incr.index <<- index*16+1:index          #incidence, previously treated
  d1.index <<- index*17+1:index            #diagnosed primary
  d2.index <<- index*18+1:index            #diagnosed secondary
  d3.index <<- index*19+1:index            #diagnosed early latent
  nsa.index <<- index*20+1:index           #not sexually active population
  d4.index <<- index*21+1:index            #diagnosed latent latent
  dr.index <<- index*22+1:index            #diagnosed primary, secondary, and early latent, previously treated
  
  pop1 <<- c(1:4,21:24) #subpop1
  pop2 <<- c(5:8,25:28) # subpop2
  pop3 <<- c(9:12,29:32) #subpop3
  pop4 <<-c(13:16,33:36 ) #subpop4
  pop5 <<- c(17:20,37:40) #subpop5
  males <<-1:20 # males
  females <<-21:40 #females
  m1<<-1:4 #M subpop1
  m2<<-5:8 #M subpop2
  m3<<-9:12 #M subpop3
  m4<<-13:16 #M subpop4
  m5<<-17:20 #M subpop5
  msw <<-1:12 #M heterosexual
  f1<<-21:24 #F subpop1
  f2<<-25:28 #F subpop2
  f3<<-29:32 #F subpop3
  f4<<-33:36 #F subpop4 (empty)
  f5<<-37:40 #F subpop5 (empty)
  y.m<<-c(1:2,5:6,9:10,13:14,17:18) #youngest age cat M
  y.m.msw<<- c(1:2, 5:6, 9:10) #yougest age cat MSW only
  y.f<<-c(21:22,25:26,29:30) #youngest age cat F
  o.m<<-c(3:4,7:8, 11:12,15:16,19:20) #oldest age cat M
  o.m.msw <<- c(3:4, 7:8, 11:12) #oldest age cat MSW only
  o.f<<-c(23:24,27:28,31:32) #oldest age cat F
  y.m.1 <<- 1:2 #youngest age cat, subpop1
  o.m.1 <<- 3:4 #old age cat, subpop1
  y.m.2 <<- 5:6 #young age cat, subpop 2
  o.m.2 <<- 7:8 # old age cat, subpop 2
  y.m.3 <<- 9:10 # young age cat, subpop 3
  o.m.3 <<- 11:12 # old age cat, subpop 3
  y.m.4 <<- 13:14 # young age cat, subpop 4
  o.m.4 <<- 15:16 # old age cat, subpop 4
  y.m.5 <<- 17:18 #young age cat, subpop 5
  o.m.5 <<- 19:20 #old age cat, subpop 5
  y.f.1 <<- 21:22 #youngest age cat, subpop1
  o.f.1 <<- 23:24 #old age cat, subpop1
  y.f.2 <<- 25:26 #young age cat, subpop 2
  o.f.2 <<- 27:28 # old age cat, subpop 2
  y.f.3 <<- 29:30 # young age cat, subpop 3
  o.f.3 <<- 31:32 # old age cat, subpop 3
  y.f.23 <<- c(25:26,29:30) # youngest age cat, subpops 2&3
  n.s.a<<-(sapply(1:((length(s.index)-4)/2),function(x){sum(n.sa[1:2+(x-1)*2])})) # pop sizes for given age, sex, subpop (sexually active only)
  n.ns.a<<-(sapply(1:((length(s.index)-4)/2),function(x){sum(n.i[1:2+(x-1)*2])})) # pop sizes for given age, sex, subpop (all)
  p.m.sa <<- c(sum(n.sa[c(m2,m3)])/sum(n.sa[c(m2,m3,m4,m5)]),sum(n.sa[c(m1,m3)])/sum(n.sa[c(m1,m3,m4,m5)]),sum(n.sa[c(m1,m2)])/sum(n.sa[c(m1,m2,m4,m5)]),1,1) #used in mixing matrix to calculate F contacts with other non-MSM subpops
  p.m.hivneg.sa <<- c(sum(n.sa[m4])/sum(n.sa[c(m2,m3,m4,m5)]),sum(n.sa[m4])/sum(n.sa[c(m1,m3,m4,m5)]),sum(n.sa[m4])/sum(n.sa[c(m1,m2,m4,m5)]),1,1)
  p.m.hivpos.sa <<- c(sum(n.sa[m5])/sum(n.sa[c(m2,m3,m4,m5)]),sum(n.sa[m5])/sum(n.sa[c(m1,m3,m4,m5)]),sum(n.sa[m5])/sum(n.sa[c(m1,m2,m4,m5)]),1,1)
  
  }

