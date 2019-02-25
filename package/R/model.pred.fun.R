###############################################################################################
### function to update estimated model parameters, run the model, and print desired outputs ###
###############################################################################################

prediction.epi <- function(theta) { # fit to infection data
  ### return log/logit transformed parameters back to untransformed values and put in required form for use in transmission model ###
  epsilon <- ilogit(c(theta["logit.epsilon.1"], theta["logit.epsilon.2"], theta["logit.epsilon.3"], theta["logit.epsilon.4"], theta["logit.epsilon.5"])) #mixing between high and low activity groups (black, other, Hispanic, MSM HIV-, MSM HIV+)
  c.min.param <- exp(c(theta["log.c.min.m.1"], theta["log.c.min.m.2"], theta["log.c.min.f.1"], theta["log.c.min.f.2"], theta["log.c.min.msm.1"], theta["log.c.min.msm.2"]))
  c.min <- c.min.fun(c.min.param) #produces a vector of minimum partner change rates 
  ctrl.pt<-update.ctrl(theta)  #internal control points of bezier curves
  bez.bc <-ctrl.pt[[2]] #get internal control points of bezier curves
  screen.bezier <- matrix(c(rep(c(ilogit(c(theta["logit.screen.m1.a"],theta["logit.screen.m1.d"])), bez.bc["m1","b"],bez.bc["m1","c"], ilogit(c(theta["logit.screen.m1.a"],theta["logit.screen.m1.d"])), bez.bc["m1","b"],bez.bc["m1","c"]),3),
                            rep(c(ilogit(c(theta["logit.screen.msm1.a"],theta["logit.screen.msm1.d"])), bez.bc["msm1","b"],bez.bc["msm1","c"], ilogit(c(theta["logit.screen.msm1.a"],theta["logit.screen.msm1.d"])), bez.bc["msm1","b"],bez.bc["msm1","c"]),2),
                            rep(c(ilogit(c(theta["logit.screen.f1.a"],theta["logit.screen.f1.d"])), bez.bc["f1","b"],bez.bc["f1","c"], ilogit(c(theta["logit.screen.f1.a"],theta["logit.screen.f1.d"])), bez.bc["f1","b"],bez.bc["f1","c"]),3),
                            rep(0,16)), ncol=4, byrow=TRUE)
  #browser()
  screen.trend <- apply(screen.bezier,1, bezier.fun) #get base annual screening rates over the calibration period
  rr.screen.o <- unname(exp(c(rep(c(log(1), theta["log.rr.screen.o.m"]),times=3), rep(c(log(1), theta["log.rr.screen.o.msm"]),times=2), rep(c(log(1), theta["log.rr.screen.o.f"]), times=5)))) #screening rr by older age/sex
  rr.screen.s<- unname(exp(rep(c(theta["log.rr.screen.m1"], log(1), theta["log.rr.screen.m3"],log(1), theta["log.rr.screen.msmhiv"], theta["log.rr.screen.f1"], log(1), theta["log.rr.screen.f3"],log(1),log(1)), each=2 ))) #screening rr by race/ethnicity/hiv status
  rr.screen.s <-rr.screen.s * rr.screen.o #screening rr by age and race/ethnicity
  rr.screen <- unname(exp(theta["log.rr.screen.ac"])) #screening rr in higher sexual activity group
  screen <- unname(t(apply(t(apply(screen.trend, 1, function(x) x*rr.screen.s)), 1, FUN=screen.fun, rr.screen=rr.screen))) #calculate actual screening rates for each group by multiplying base rates by rr
  screen[screen<0.0001]<-0  # make sure that screening rate isn't <0
  screen <-rbind(screen[rep(1, times=(cal.start-11)),], screen, screen[rep((nrow(screen)), times=10),]) ##generate screening rates for the entire model period (assume screening rates pre-calibration start=rate at start of calibration); time trend starts 10 years before cal start
  b <- ilogit(c(theta["logit.b.m"], theta["logit.b.f"], theta["logit.b.msm"])) #transmission rates
  delta <- 365/((exp(theta["log.dur.incub"]))) #1/incubation period
  dur.inf <- c(exp(theta["log.dur.prim"]), exp(theta["log.dur.sec"]),365-(exp(theta["log.dur.prim"])+exp(theta["log.dur.sec"]))) #infectious durations from primary, secondary, and early latent sypilis (early latent dur = 365 d - dur prim - dur sec)
  gamma <- c(365/(dur.inf[1]), 365/(dur.inf[2]), 365/(dur.inf[3]) ) #1/dur infectious
  dur.imm <- exp(c(theta["log.dur.imm.inf"], (theta["log.dur.imm.inf"]+theta["log.dur.imm.early"]), theta["log.dur.immune"])) #duration of protective immunity after treatment for p&s, early latent, or late latent syphilis
  p.trt.1 <- ilogit(c(theta["logit.p.trt.prim.m"], theta["logit.p.trt.prim.msm"], theta["logit.p.trt.prim.msmhiv"], theta["logit.p.trt.prim.f"])) #treatment rate for primary syphilis (male, hiv-msm, hiv+msm, female)
  p.trt.2 <- ilogit(c(theta["logit.p.trt.sec.m"], theta["logit.p.trt.sec.msm"], theta["logit.p.trt.sec.msmhiv"], theta["logit.p.trt.sec.f"]))  #treatment rate for seconary syphilis (male, hiv-msm, hiv+msm, female)
  p.trt.3 <- ilogit(c(theta["logit.p.trt.lat.m"], theta["logit.p.trt.lat.msm"], theta["logit.p.trt.lat.msmhiv"], theta["logit.p.trt.lat.f"]))  #treatment rate for early latent syphilis (male, hiv-msm, hiv+msm, female)
  p.trt.4 <- ilogit(c(theta["logit.p.trt.late.m"], theta["logit.p.trt.late.msm"], theta["logit.p.trt.late.msmhiv"], theta["logit.p.trt.late.f"]))  #treatment rate for late latent syphilis (male, hiv-msm, hiv+msm, female)
  theta.param <- ilogit(c(theta["logit.theta.1"], theta["logit.theta.2"], theta["logit.theta.3"], theta["logit.theta.4"], theta["logit.theta.5"], theta["logit.theta.6"], theta["logit.theta.7"],theta["logit.theta.8"], rep(logit(1),2))) #subpopulation assortativity
  theta.param <- matrix(theta.param,ncol=2) #col1=males, col2=females
  theta.hiv <- ilogit(theta["logit.theta.hiv"]) #proportion of seroconcordant partnerships,for HIV+ MSM
  pi.all<- ilogit(c(rep(c(theta["logit.pi.m"],theta["logit.pi.f"]),3),rep(theta["logit.pi.msm"],4))) #age assortativity
  rr.rep.symp.m <- ilogit(theta["logit.rr.rep.symp.m"]) #relative risk case is reported if seeks treatment and male
  rr.rep.symp.f <-ilogit(theta["logit.rr.rep.symp.f"]) #realtive risk case is reported if seeks treatment and female
  rep.bc<-ctrl.pt[[1]] #get internal control points for reporting bezier curve
  rep.bezier <- matrix(rep(c(ilogit(theta["logit.rep.a"]), ilogit(theta["logit.rep.d"]),rep.bc[1],rep.bc[2]),40), ncol=4, byrow=TRUE) #get 4 bezier points for calc reporting prob in each subgroup
  rep.trend <- apply(rep.bezier,1, bezier.fun) #baseline reporting probability over calibration period
  rep.trend[rep.trend>1]<-1  # make sure that reporting rate isn't >1
  rep <-rbind(rep.trend[rep(1, times=(cal.start-11)),], rep.trend, rep.trend[rep((nrow(rep.trend)), times=10),])  #expand reporting prob to cover entire model run time, assume reporting pre-calibration period = reporting at start of cal; time trend starts 10 years prior to cal start
  rep.symp<-matrix(c(rep[,males]*rr.rep.symp.m, rep[,females]*rr.rep.symp.f), ncol=40) #reporting rate if male or female case who actively seeks treatment
  rep<-rep[,1] #reporting rate if case identified by screening
  behav.trend<-behav.fun(ilogit(theta["logit.behav.lin"])) # behav time trend representing changing condom use/behaviour in MSM
  behav <-c(rep(1, times=(cal.start-11)), behav.trend, behav.trend[rep(length(behav.trend), times=10)]) #expand to cover entire model run time, assume behav pre-calibration period = value at start of calibration; trend starts 10 years prior to cal start
  rp.all <- exp(c(theta["log.rp.1.1.1.1"],theta["log.rp.1.1.2.1"],theta["log.rp.1.2.1.1"],theta["log.rp.1.2.2.1"], #relative rates of partner change in different population groups
                  theta["log.rp.1.1.1.2"],theta["log.rp.1.1.2.2"],theta["log.rp.1.2.1.2"], theta["log.rp.1.2.2.2"],
                  theta["log.rp.2.1.1.1"],theta["log.rp.2.1.2.1"],theta["log.rp.2.2.1.1"],theta["log.rp.2.2.2.1"],
                  theta["log.rp.2.1.1.2"],theta["log.rp.2.1.2.2"],theta["log.rp.2.2.1.2"], theta["log.rp.2.2.2.2"],
                  log(1),theta["log.rp.3.1.2.1"],log(1),theta["log.rp.3.2.2.1"],
                  log(1),theta["log.rp.3.1.2.2"],log(1), theta["log.rp.3.2.2.2"],
                  log(1),theta["log.rp.4.1.2.1"],log(1),log(1),
                  log(1),theta["log.rp.4.1.2.2"],log(1), log(1),
                  theta["log.rp.5.1.1.1"],theta["log.rp.5.1.2.1"],log(1),log(1), ##rp for HIV+ MSM - will need to update##
                  theta["log.rp.5.1.1.2"],theta["log.rp.5.1.2.2"],log(1), log(1)
  ))
  pred <- mixing(epsilon, pi.all, theta.param, c.min, rp.all, theta.hiv) #calculate mixing matrix based on current parameter estimates
  age.dist.all<-prop.table(apply(part.all, c(3:5), sum), c(1,3)) #get proportion of partnerships with same age group, by sex
  age.dist.all<-c(diag(age.dist.all[,,1]),diag(age.dist.all[,,2])) #proportion of partners of same age, entire pop: M-age1, age2, F-age1, age2
  age.dist.sub<-prop.table(apply(part.all, c(3:6), sum),c(1,3,4) )  #get proportion of partnerships with same age group, by subpopulation
  d.m=c()  
  d.f=c()
  for (i in 1:5){
    d.m[(2*i-1):(2*i)]<-diag(age.dist.sub[,,1,i])
    d.f[(2*i-1):(2*i)]<-diag(age.dist.sub[,,2,i])
  }
  age.dist.all<-c(age.dist.all, d.m[1:(3*2)], d.f[1:(3*2)],d.m[7:10])  #add in estimate for MSM 
  s.dist.m <- part.all.m[1:3,1:3]/apply(part.all.m[1:3,1:3],1,sum)  #proportion of partners of same or other subpopulation, M (excluding MSM)
  s.dist.f <- part.all.f[1:3,1:3]/apply(part.all.f[1:3,1:3],1,sum) #proportion of parters of same or other subpopulation, F
  s.dist.msm <- part.all.m[4:5,4:5]/apply(part.all.m[4:5,4:5],1,sum)  #proportion of MSM partnerships with same HIV status (for HIV- and HIV+)
  pred.s.dist <- c(diag(s.dist.m),diag(s.dist.f),diag(s.dist.msm)) #partners of same subpopulation, M and F
  replacement <- matrix(0L, nrow = dim(screen[102:105,])[1], ncol = dim(screen[102:105,])[2]) # turn off screening in 2013 and afterwards
  screen[102:105,] <- replacement
  params <-list(b=b,delta=delta, gamma=gamma,p.trt.1=p.trt.1, p.trt.2=p.trt.2, p.trt.3=p.trt.3,p.trt.4=p.trt.4, rep=rep,rep.symp=rep.symp, screen=screen, behav=behav,dur.imm=dur.imm, ct.data.years=ct.data.years, p.ct.primsec=p.ct.primsec, p.ct.el = p.ct.el, p.s.1=p.s.1, p.s.2 = p.s.2, p.s.3 = p.s.3, d=d, fileName="tintoTrt.txt", fileNameET="texitRateTrt.txt", fileNameES="texitRateScr.txt", fileNameA="talpha.txt", fileNameC="tcovps.txt") #parameters used by transmission model
  ### run the transmission model ###
  out.cpp<-syphSim(params,
                   tstep,
                   cm.list,
                   p.abx,
                   rep.count, 
                   model.end+1,
                   yinit, 
                   n.sa,
                   births,
                   births.sa, 
                   births.nsa,
                   aging,
                   aging.nsa)
  sol<-out.cpp$out
  ct <- out.cpp$ct
  
  prev<-model.prev(sol)  #calculate desired outputs for calibration, function stored in model.out.fun.R
  prev.cf <- NULL #initialize in case we're not running the counterfactual
  if(showCounterfactual == TRUE) { #counterfactual for no contact tracing, based on imputed coverage
    screen.cf <- screen
    screen.cf[(model.end-ct.data.years):(model.end-1),] <- screen[(model.end-ct.data.years):(model.end-1),] - ct*(1/d)
      #screen[100:105,] - ct * screen[100:105,]
    screen.cf[model.end:nrow(screen.cf),] <- screen.cf[(model.end-1),]
    replacement <- matrix(0L, nrow = dim(screen.cf[100:105,])[1], ncol = dim(screen.cf[100:105,])[2])
    #replacement <- matrix(0L, nrow = dim(screen.cf)[1], ncol = dim(screen.cf)[2])
    screen.cf[100:105,] <- replacement
    #browser()
    screen.cf[is.nan(screen.cf)] <- 0
    # screen.cf <- matrix(1, nrow=nrow(screen.cf), ncol=ncol(screen.cf))
    params.cf <-
      list(
        b = b,
        delta = delta,
        gamma = gamma,
        p.trt.1 = p.trt.1,
        p.trt.2 = p.trt.2,
        p.trt.3 = p.trt.3,
        p.trt.4 = p.trt.4,
        rep = rep,
        rep.symp = rep.symp,
        screen = screen.cf,
        behav = behav,
        dur.imm = dur.imm,
        ct.data.years = ct.data.years,
        p.ct.primsec = p.ct.primsec,
        p.ct.el = p.ct.el,
        p.s.1 = p.s.1,
        p.s.2 = p.s.2,
        p.s.3 = p.s.3,
        d = d,
        fileName="tintoTrtCF.txt", fileNameET="texitRateTrtCF.txt", fileNameES="texitRateScrCF.txt", fileNameA="talphaCF.txt", fileNameC="tcovpsCF.txt"
      ) #parameters used by transmission model
    ### run the transmission model ###
    out.cpp.cf<-syphSim(params.cf,
                     tstep,
                     cm.list,
                     p.abx,
                     rep.count, 
                     model.end+1,
                     yinit, 
                     n.sa,
                     births,
                     births.sa, 
                     births.nsa,
                     aging,
                     aging.nsa)
    sol.cf<-out.cpp.cf$out

    prev.cf<-model.prev(sol.cf)
  }
  #browser()
  list(subpop.assort=pred.s.dist, age.assort=age.dist.all, prev=prev, prev.cf=prev.cf)
}



