#########################################################
### functions for calculating time-varying parameters ###
#########################################################
# bezier function for screening rate/reporting, constrained to be between 0 and 1
# function takes 4 parameters: 1 = start point (value at t=0), 2=end point (value at t=end), 3&4=internal control points
# produces annual estimates of screening/reporting rate for length specified (the calibration period)

bezier.fun<- function(screen.bezier){
  bez.a<-screen.bezier[1]
  bez.d<-screen.bezier[2]
  bez.b <- screen.bezier[3]
  bez.c <- screen.bezier[4]
  y0=bez.a
  y1=bez.b
  y2=bez.c
  y3=bez.d
  p0 =  y0
  p1 = ( -5*y0 + 18*y1 -  9*y2 + 2*y3) / 6
  p2 = (  2*y0 -  9*y1 + 18*y2 - 5*y3) / 6
  p3 = y3
  p<-c(p0,p1,p2,p3)
  t<-seq(0,1, length.out=(cal.period+10)) #let time trend start 10 years before calibration start
  #browser()
  bezier::bezier(t,p)  #this is the function from the bezier package - Hmisc also has a bezier function that is NOT the one to use
}


#calculate control points b and c for bezier curves describing screening and reporting
#see model technical appendix for additonal details 
update.ctrl <- function(theta) {
  #internal control points for screening rate in females
  screen.f1.b = ilogit(theta["logit.screen.f1.a"]) + (ilogit(theta["logit.screen.f1.d"])-ilogit(theta["logit.screen.f1.a"]))*ilogit(theta["logit.rand.screen.f1.b"])
  screen.f1.c = screen.f1.b + (ilogit(theta["logit.screen.f1.d"])-screen.f1.b)*ilogit(theta["logit.rand.screen.f1.c"])
  #internal control points for screening rate in males
  screen.m1.b = ilogit(theta["logit.screen.m1.a"]) + (ilogit(theta["logit.screen.m1.d"])-ilogit(theta["logit.screen.m1.a"]))*ilogit(theta["logit.rand.screen.m1.b"])
  screen.m1.c = screen.m1.b + (ilogit(theta["logit.screen.m1.d"])-screen.m1.b)*ilogit(theta["logit.rand.screen.m1.c"])
  #internal control points for screening rate in MSM
  screen.msm1.b = ilogit(theta["logit.screen.msm1.a"]) + (ilogit(theta["logit.screen.msm1.d"])-ilogit(theta["logit.screen.msm1.a"]))*ilogit(theta["logit.rand.screen.msm1.b"])
  screen.msm1.c = screen.msm1.b + (ilogit(theta["logit.screen.msm1.d"])-screen.msm1.b)*ilogit(theta["logit.rand.screen.msm1.c"])
  screen.bez.ctrl <-matrix(c(screen.m1.b, screen.m1.c,
                             screen.msm1.b, screen.msm1.c,
                             screen.f1.b, screen.f1.c), ncol=2, byrow=T)
  rownames(screen.bez.ctrl)<-c("m1","msm1","f1")
  colnames(screen.bez.ctrl)<-c("b","c")
  #internal control points for reporting parameter
  rep.b = ilogit(theta["logit.rep.a"]) + (ilogit(theta["logit.rep.d"])-ilogit(theta["logit.rep.a"]))*ilogit(theta["logit.rand.rep.b"])
  rep.c = rep.b + (ilogit(theta["logit.rep.d"])-rep.b)*ilogit(theta["logit.rand.rep.c"])
  rep.bez.bc <<- c(rep.b, rep.c)     
  
  return(list(rep.bez.bc,screen.bez.ctrl))
}

#calcuate internal control points for bezier curves for plotting model priors
prior.ctrl <- function(bez) {   #update control points b and c for bezier curves
  bez.b <- bez[1] + (bez[2]-bez[1])*bez[3]
  bez.c <- bez.b + (bez[2]-bez.b)*bez[4]
  prior <- unlist(c(bez[1], bez[2], b=bez.b, c=bez.c) )
  return(prior)
}

#calculate screening rate from base screening rate x rr screening in high sexual activity group, for each row of matrix
screen.fun <-function(screen.param,rr.screen) {
  x=NULL
  for (i in 1:length(screen.param)) {
    x<-c(x,screen.param[i], screen.param[i]*rr.screen)
  }
  return(x)
}

#calculate time-varying tranmission relative risk for MSM
behav.fun<- function(behav.m){
  # t<-seq(0,cal.period+9, by=1) #let time trend start 10 years before model calibration start
  t <- seq(0,1,length.out = cal.period+9)
  return(behav.m*t)
}
