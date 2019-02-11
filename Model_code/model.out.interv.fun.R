################################################################################
### function to calculate outputs used for comparing impact of interventions ###
################################################################################

model.prev <- function(sol, interv) {  #sol is the output from the transmissio model, a matrix of model compartment sizes at each time step
  #calculate total population size, sexually active (sa) and total
  pop.size.sa <- sol[,1+s.index] + sol[,1+e.index] + sol[,1+prim.index] + sol[,1+sec.index] +  sol[,1+early.index] +
     sol[,1+latent.index] +sol[,1+treated.inf.index] + sol[,1+treated.early.index] + sol[,1+treated.late.index] + sol[,1+sr.index] + 
     sol[,1+er.index] + sol[,1+primr.index] + sol[,1+secr.index] +  sol[,1+earlyr.index] + sol[,1+latentr.index] 
  pop.size.all <- pop.size.sa + sol[,1+nsa.index]
  #calculate prevalent infections
  inf.ps <- sol[,1+prim.index]+sol[,1+sec.index] + sol[,1+primr.index]+sol[,1+secr.index]  #primary and secondary infections
  inf.latent <- sol[,1+early.index] + sol[,1+earlyr.index] #early latent infections
  inf.late <- sol[,1+latent.index] + sol[,1+latentr.index] #late latent infection
  inf.total <- inf.ps+inf.latent #all early latent infections
  p.early.rep <- (sol[,1+primr.index] + sol[,1+secr.index] + sol[,1+earlyr.index]) / inf.total
  pop.subpop.age <- colSums(matrix(tail(pop.size.all,1), nrow=2))
  pop.sub.adj <- rep(pop.subpop.age[11:16],2) #pop size by sex, subpop, redistributing MSM for calculating reported rates
  #prevalence of early syphilis
  prev.total <- cbind(rowSums(inf.total[,m1])/sum(n.i[m1]), rowSums(inf.total[,m2])/sum(n.i[m2]),rowSums(inf.total[,m3])/sum(n.i[m3]),rowSums(inf.total[,m4])/sum(n.i[m4]),  
                      rowSums(inf.total[,m5])/sum(n.i[m5]),rowSums(inf.total[,f1])/sum(n.i[f1]),rowSums(inf.total[,f2])/sum(n.i[f2]),rowSums(inf.total[,f3])/sum(n.i[f3]))
  prev.tot.sex <- cbind(rowSums(inf.total[,males])/sum(n.i[males]), rowSums(inf.total[,females])/sum(n.i[females]))
  
  #late latent syphilis prevalence
  prev.late <- cbind(rowSums(inf.late[,m1])/sum(n.i[m1]), rowSums(inf.late[,m2])/sum(n.i[m2]),rowSums(inf.late[,m3])/sum(n.i[m3]),rowSums(inf.late[,m4])/sum(n.i[m4]),
                     rowSums(inf.late[,m5])/sum(n.i[m5]),rowSums(inf.late[,f1])/sum(n.i[f1]),rowSums(inf.late[,f2])/sum(n.i[f2]),rowSums(inf.late[,f3])/sum(n.i[f3]))
  prev.late.sex <- cbind(rowSums(inf.late[,males])/sum(n.i[males]), rowSums(inf.late[,females])/sum(n.i[females]))
  ###calculate incidence for period covered by calibration data
  inc <- annual.rates(sol[,1+inc.index]) #annual.rates assumes model is outputting ANNUAL values -- would need to change function if producing weekly/monthly outputs
  inc.s.a <- (sapply(1:((ncol(inc)-8)/2),function(x){rowSums(inc[,1:2+(x-1)*2])})) #incident CASES for given age, sex, and subpopulation group
  n.inc <<- as.vector(inc.s.a[(cal.start+1):(cal.start+cal.period),] ) #incident cases for calibration period
  inc.age.sex <- cbind(rowSums(inc[,y.m])/sum(n.i[y.m]), rowSums(inc[,o.m])/sum(n.i[o.m]), rowSums(inc[,y.f])/sum(n.i[y.f]), rowSums(inc[,o.f])/sum(n.i[o.f]), rowSums(inc[,c(y.m.4,y.m.5)])/sum(n.i[c(y.m.4, y.m.5)]))
  cum.inc.m <- sum(inc[(interv.start):(interv.end),males]) #incidence in males 
  cum.inc.msw <- sum(inc[(interv.start):(interv.end),msw])#overall incidence in msw
  cum.inc.msm <- sum(inc[(interv.start):(interv.end),c(m4, m5)]) #incidence in msm 
  cum.inc.f <- sum(inc[(interv.start):(interv.end),females]) #overall incidence in females
  cum.inc.tot <- cum.inc.m + cum.inc.f
  #if(interv==0) {
    #default.cum.inc.m <<- cum.inc.m
    #default.cum.inc.msw <<- cum.inc.msw
    #default.cum.inc.msm <<- cum.inc.msm
    #default.cum.inc.f <<- cum.inc.f
    #default.cum.inc.tot <<- cum.inc.tot
  #} else {
    
  #}
  #browser() 
  diag.p <- annual.rates(sol[,1+d1.index]) #reported primary cases
  diag.s <- annual.rates(sol[,1+d2.index]) #reported secondary cases
  diag.el <- annual.rates(sol[,1+d3.index]) #reported early latent cases
  diag.ll <- annual.rates(sol[,1+d4.index]) #reported late latent cases
  diag.early <- diag.p + diag.s + diag.el #reported early syphilis cases
  diag.early.s.a <- (sapply(1:((ncol(diag.early)-8)/2),function(x){rowSums(diag.early[,1:2+(x-1)*2])})) #reported early CASES for given age, sex, and subpopulation group
  diag.early.s.a.adj <-diag.early.s.a[,1:6] + cbind(diag.early.s.a[,7:8]*p.s.1, diag.early.s.a[,7:8]*p.s.2, diag.early.s.a[,7:8]*p.s.3) + cbind(diag.early.s.a[,9:10]*p.s.1, diag.early.s.a[,9:10]*p.s.2, diag.early.s.a[,9:10]*p.s.3)#redistribute MSM cases among other male subpops
  diag.early.adj <- cbind(diag.early.s.a.adj,diag.early.s.a[,11:16]) #append female cases to male cases with MSM distributed among other male cases
  diag.early.adj.rate <-100000*t(t(diag.early.adj)/pop.sub.adj) #reported cases rates adjusted for MSM
  
  cal.start.diag <- cal.start+cal.period - nrow(diag.age.sex.dat) + 1  #calculate years of reported case data (by age and sex) based on input (all data end in year 2016)
  diag.age.sex.fit.m <- 100000*c(rowSums(diag.early[((cal.start.diag):(cal.start+cal.period)),y.m])/sum(n.i[y.m]),  #reported early syphilis cases, young and old males
                                rowSums(diag.early[((cal.start.diag):(cal.start+cal.period)),o.m])/sum(n.i[o.m]))  
  diag.age.sex.fit.f.y <- 100000*(rowSums(diag.early[((cal.start.diag):(cal.start+cal.period)),y.f])/sum(n.i[y.f])) #reported early syphilis cases, young females
  diag.age.sex.fit.f.o <- 100000*(rowSums(diag.early[((cal.start.diag):(cal.start+cal.period)),o.f])/sum(n.i[o.f])) #reported early syphilis cases, old females
  
  diag.age.sex.fit.ll <- 100000*c(rowSums(diag.ll[((cal.start.diag):(cal.start+cal.period)),y.m])/sum(n.i[y.m]), #reported late syphilis cases, by sex and age group
                                 rowSums(diag.ll[((cal.start.diag):(cal.start+cal.period)),o.m])/sum(n.i[o.m]), 
                                 rowSums(diag.ll[((cal.start.diag):(cal.start+cal.period)),y.f])/sum(n.i[y.f]),
                                 rowSums(diag.ll[((cal.start.diag):(cal.start+cal.period)),o.f])/sum(n.i[o.f]))  
### calculate reported case relative rate by race/ethnicity, using overall rates as comparator
### using average of last 5 years for model fitting due to instability in the estimates from year to year
  cal.start.rr <- cal.start+cal.period+1 -5 #using average of last 5 years for calibration
  diag.age.sex.rate <- 100000*c(rowSums(diag.early[((cal.start.rr):(cal.start+cal.period)),y.m])/sum(n.i[y.m]), #overall rates by age/sex for denominator
                                rowSums(diag.early[((cal.start.rr):(cal.start+cal.period)),o.m])/sum(n.i[o.m]),  
                                rowSums(diag.early[((cal.start.rr):(cal.start+cal.period)),y.f])/sum(n.i[y.f]), 
                                rowSums(diag.early[((cal.start.rr):(cal.start+cal.period)),o.f])/sum(n.i[o.f])) 
  rr.diag.y.m.1 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),1])/ mean(diag.age.sex.rate[1:5])  ## relative diagnoses, average 2012-2016
  rr.diag.o.m.1 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),2])/ mean(diag.age.sex.rate[6:10])
  rr.diag.y.m.3 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),5])/ mean(diag.age.sex.rate[1:5])
  rr.diag.o.m.3 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),6])/ mean(diag.age.sex.rate[6:10])
  rr.diag.y.f.1 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),7])/ mean(diag.age.sex.rate[11:15])
  rr.diag.o.f.1 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),8])/ mean(diag.age.sex.rate[16:20])
  rr.diag.y.f.3 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),11])/ mean(diag.age.sex.rate[11:15])
  rr.diag.o.f.3 <- mean(diag.early.adj.rate[(cal.start.rr):(cal.start+cal.period),12])/ mean(diag.age.sex.rate[16:20])
  fit.diag.rr <- c(rr.diag.y.m.1, rr.diag.o.m.1, rr.diag.y.m.3, rr.diag.o.m.3, rr.diag.y.f.1, rr.diag.o.f.1, rr.diag.y.f.3, rr.diag.o.f.3)
  
### calculate proportion of cases in MSM
  cal.start.msm <- cal.start+ cal.period +1 - nrow(msm.dat[,c("pMSM_y", "pMSM_o")])  #calculate number of years of data available
  p.diag.msm <- c(rowSums(diag.early[(cal.start.msm):(cal.start+cal.period), c(y.m.4, y.m.5)])/ rowSums(diag.early[(cal.start.msm):(cal.start+cal.period), y.m]),
                  rowSums(diag.early[(cal.start.msm):(cal.start+cal.period), c(o.m.4, o.m.5)])/ rowSums(diag.early[(cal.start.msm):(cal.start+cal.period), o.m])) #proportion of total MALE diagnosed cases in MSM (2007-2016)
### calculate proportion of MSM cases with HIV coinfection  
  cal.start.hiv <- cal.start + cal.period + 1 - nrow(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"])))) #years of data for pMSM with HIV
  p.diag.hiv <- c(rowSums(diag.early[(cal.start.hiv):(cal.start+cal.period), y.m.5])/ rowSums(diag.early[(cal.start.hiv):(cal.start+cal.period), c(y.m.4,y.m.5)]),
                  rowSums(diag.early[(cal.start.hiv):(cal.start+cal.period), o.m.5])/ rowSums(diag.early[(cal.start.hiv):(cal.start+cal.period), c(o.m.4,o.m.5)])) #proportion of MSM cases in HIV+ (2007-2016)
#calculate distribution of early syphilis cases by stage at diagnosis
  #proportion of early cases diagnosed with secondary syphilis
  p.diag.sec <- c(rowSums(diag.s[(cal.start.diag):(cal.start+cal.period), y.m]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), y.m]),
                  rowSums(diag.s[(cal.start.diag):(cal.start+cal.period), o.m]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), o.m]),
                  rowSums(diag.s[(cal.start.diag):(cal.start+cal.period), y.f]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), y.f])
                  )
 #proportion of early cases diagnosed with early latent syphilis 
  p.diag.el <- c(rowSums(diag.el[(cal.start.diag):(cal.start+cal.period), y.m]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), y.m]), 
                 rowSums(diag.el[(cal.start.diag):(cal.start+cal.period), o.m]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), o.m]),
                 rowSums(diag.el[(cal.start.diag):(cal.start+cal.period), y.f]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), y.f])
                  )
  #proportion of all diagnoses syphilis cases diagnosed with early (primary, secondary, or early latent) syphilis
  p.diag.early <- mean(rowSums(diag.early[(cal.start.diag):(cal.start+cal.period),]) / (rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), ]) + rowSums(diag.ll[(cal.start.diag):(cal.start+cal.period), ])))
  
  #proportion of diagnosed early syphilis cases in those with prior infection
  diag.early.r <- annual.rates(sol[,1+dr.index]) #reported early syphilis cases, prior infection
  p.prior.diag <-  c(rowSums(diag.early.r[(cal.start.diag):(cal.start+cal.period), males]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), males]), 
                rowSums(diag.early.r[(cal.start.diag):(cal.start+cal.period), females]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), females]),
                rowSums(diag.early.r[(cal.start.diag):(cal.start+cal.period), c(m4,m5)]) /rowSums(diag.early[(cal.start.diag):(cal.start+cal.period), c(m4,m5)])
  )
  
  #proportion of incident syphilis cases in those with prior infection
  inc.r <- annual.rates(sol[,1+incr.index]) #incident syphilis cases, prior infection
  p.prior.inc <-  c(rowSums(inc.r[(cal.start.diag):(cal.start+cal.period), males]) /rowSums(inc[(cal.start.diag):(cal.start+cal.period), males]), 
                     rowSums(inc.r[(cal.start.diag):(cal.start+cal.period), females]) /rowSums(inc[(cal.start.diag):(cal.start+cal.period), females]),
                     rowSums(inc.r[(cal.start.diag):(cal.start+cal.period), c(m4,m5)]) /rowSums(inc[(cal.start.diag):(cal.start+cal.period), c(m4,m5)])
  )
  
  list(early.inf.rate.m= diag.age.sex.fit.m, early.inf.rate.f.y= diag.age.sex.fit.f.y, early.inf.rate.f.o = diag.age.sex.fit.f.o, diag.rr = fit.diag.rr, p.diag.msm = p.diag.msm, n.inc=n.inc, p.diag.hiv=p.diag.hiv, p.diag.sec=p.diag.sec, p.diag.el=p.diag.el, diag.late=diag.age.sex.fit.ll, p.diag.early=p.diag.early, 
       p.prior.inc = p.prior.inc, p.prior.diag=p.prior.diag, cum.inc.m=cum.inc.m, cum.inc.f=cum.inc.f, cum.inc.msm=cum.inc.msm, cum.inc.tot=cum.inc.tot)
}

