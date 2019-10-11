###########################################################################
### visualize prior and posterior distributions and calibration targets ###
###########################################################################

#run this after calibration, with burned, trimmed, and merged trace

library(gridExtra)
library(grid)
library(Hmisc)
library(R.utils)

bezier <- bezier::bezier 

get_legend<-function(myggplot){   # function to get legend from plot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blankPlot <- ggplot()+geom_blank(aes(1,1)) +  # make a blank plot to use as a placeholder where necessary
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )


# plot calibration targets and priors/posteriors 

plot.posteriors <- function(post.sample, output_dir) {
  #cal.period <- st_env$cal.period
  #n.i <- st_env$n.i
  #y.m.1 <- st_env$y.m.1
  #y.m.2 <- st_env$y.m.2
  #y.m.3 <- st_env$y.m.3
  #y.m.4 <- st_env$y.m.4
  #y.m.5 <- st_env$y.m.5
  #o.m.1 <- st_env$o.m.1
  #o.m.2 <- st_env$o.m.2
  #o.m.3 <- st_env$o.m.3
  #o.m.4 <- st_env$o.m.4
  #o.m.5 <- st_env$o.m.5
  #y.f.1 <- st_env$y.f.1
  #y.f.2 <- st_env$y.f.2
  #y.f.3 <- st_env$y.f.3
  #o.f.1 <- st_env$o.f.1
  #o.f.2 <- st_env$o.f.2
  #o.f.3 <- st_env$o.f.3
  #end.year <- st_env$end.year
  #start.year <- st_env$start.year
  #prep model outputs for plotting
  # here y=20-44y, o=45-64y, m=male, f=female, 1=black, 2=other, 3=Hispanic, msm=hiv- men who have sex with men, msmhiv= hiv+ msm
  # pred <- as.data.frame(post.sample$outputs)
	pred <- t(as.data.frame(apply(post.sample$outputs, 1, unlist)))
  diag.len <- nrow(diag.age.sex.dat)
  msm.len <- nrow(msm.dat) #years of data for proportion of male cases reported as msm
  hiv.len <- nrow(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"])))) #years of data for hiv coinfection
  
  out.subpop <- reshape::melt(as.matrix(subset(pred, select=subpop.assort1:subpop.assort8))) #model subpopulation assortativity 
  out.age <- reshape::melt(as.matrix(subset(pred, select=age.assort1:age.assort4))) #model age assortativity
  out.diag.y.m <- reshape::melt(as.matrix(subset(pred, select=prev.early.inf.rate.m1:(prev.early.inf.rate.m1+diag.len-1))))  #reported early syph cases by age-sex
  out.diag.o.m <- reshape::melt(as.matrix(subset(pred, select=(prev.early.inf.rate.m1+diag.len):(prev.early.inf.rate.m1+diag.len*2-1))))  
  out.diag.y.f <- reshape::melt(as.matrix(subset(pred, select=prev.early.inf.rate.f.y1:(prev.early.inf.rate.f.y1+diag.len-1)))) 
  out.diag.o.f <- reshape::melt(as.matrix(subset(pred, select=prev.early.inf.rate.f.o1:(prev.early.inf.rate.f.o1+diag.len-1)))) 
  out.diaglate.y.m <- reshape::melt(as.matrix(subset(pred, select=prev.diag.late1:(prev.diag.late1+diag.len-1))))  #reported late syph cases by age-sex
  out.diaglate.o.m <- reshape::melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len):(prev.diag.late1+diag.len*2-1))))  
  out.diaglate.y.f <- reshape::melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len*2):(prev.diag.late1+diag.len*3-1))))
  out.diaglate.o.f <- reshape::melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len*3):(prev.diag.late1+diag.len*4-1))))  
  syph.ratio.y.m <- out.diag.y.m$value/(out.diaglate.y.m$value+ out.diag.y.m$value) #proportion of early cases to all cases by age-sex
  syph.ratio.o.m <- out.diag.o.m$value/(out.diaglate.o.m$value+ out.diag.o.m$value)
  syph.ratio.y.f <- out.diag.y.f$value/(out.diaglate.y.f$value+ out.diag.y.f$value)
  syph.ratio.o.f <- out.diag.o.f$value/(out.diaglate.o.f$value+ out.diag.o.f$value)
  
  if(exists('showCounterfactual') && showCounterfactual == TRUE) {
  out.diag.y.m.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.early.inf.rate.m1:(prev.cf.early.inf.rate.m1+diag.len-1))))  #reported early syph cases by age-sex
  out.diag.o.m.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.early.inf.rate.m1+diag.len):(prev.cf.early.inf.rate.m1+diag.len*2-1))))  
  out.diag.y.f.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.early.inf.rate.f.y1:(prev.cf.early.inf.rate.f.y1+diag.len-1)))) 
  out.diag.o.f.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.early.inf.rate.f.o1:(prev.cf.early.inf.rate.f.o1+diag.len-1)))) 
  out.diaglate.y.m.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.diag.late1:(prev.cf.diag.late1+diag.len-1))))  #reported late syph cases by age-sex
  out.diaglate.o.m.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.diag.late1+diag.len):(prev.cf.diag.late1+diag.len*2-1))))  
  out.diaglate.y.f.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.diag.late1+diag.len*2):(prev.cf.diag.late1+diag.len*3-1))))
  out.diaglate.o.f.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.diag.late1+diag.len*3):(prev.cf.diag.late1+diag.len*4-1))))  
  syph.ratio.y.m.cf <- out.diag.y.m.cf$value/(out.diaglate.y.m.cf$value+ out.diag.y.m.cf$value) #proportion of early cases to all cases by age-sex
  syph.ratio.o.m.cf <- out.diag.o.m.cf$value/(out.diaglate.o.m.cf$value+ out.diag.o.m.cf$value)
  syph.ratio.y.f.cf <- out.diag.y.f.cf$value/(out.diaglate.y.f.cf$value+ out.diag.y.f.cf$value)
  syph.ratio.o.f.cf <- out.diag.o.f.cf$value/(out.diaglate.o.f.cf$value+ out.diag.o.f.cf$value)
  }
  out.rr.diag <- reshape::melt(as.matrix(subset(pred, select=prev.diag.rr1:prev.diag.rr8))) ## reported case relative risk (pooled estimates for last 5 years)
  out.p.msm.y <- reshape::melt(as.matrix(subset(pred, select=prev.p.diag.msm1:(prev.p.diag.msm1+msm.len-1)))) #proportion of male cases in young MSM
  out.p.msm.o <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.msm1+msm.len):(prev.p.diag.msm1+msm.len*2-1)))) #proportion of male cases in old MSM
  out.p.hiv.y <- reshape::melt(as.matrix(subset(pred, select=prev.p.diag.hiv1:(prev.p.diag.hiv1+hiv.len-1)))) #proportion of young MSM cases with HIV
  out.p.hiv.o <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.hiv1+hiv.len):(prev.p.diag.hiv1+hiv.len*2-1)))) #proportion of old MSM cases with HIV
  out.p.sec.y.m <- reshape::melt(as.matrix(subset(pred, select=prev.p.diag.sec1:(prev.p.diag.sec1+diag.len-1)))) #proportion of early cases that are secondary
  out.p.sec.o.m <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.sec1+diag.len):(prev.p.diag.sec1+diag.len*2-1)))) #proportion of early cases that are secondary
  out.p.sec.y.f <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.sec1+diag.len*2):(prev.p.diag.sec1+diag.len*3-1)))) #proportion of early cases that are secondary
  out.p.el.y.m <- reshape::melt(as.matrix(subset(pred, select=prev.p.diag.el1:(prev.p.diag.el1+diag.len-1)))) #proportion of early cases that are early latent
  out.p.el.o.m <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.el1+diag.len):(prev.p.diag.el1+diag.len*2-1)))) #proportion of early cases that are early latent
  out.p.el.y.f <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.el1+diag.len*2):(prev.p.diag.el1+diag.len*3-1)))) #proportion of early cases that are early latent
  out.p.early <- reshape::melt(as.matrix(subset(pred, select=(prev.p.diag.early)))) #proportion of early cases that are early latent
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  out.rr.diag.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.diag.rr1:prev.cf.diag.rr8))) ## reported case relative risk (pooled estimates for last 5 years)
  out.p.msm.y.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.p.diag.msm1:(prev.cf.p.diag.msm1+msm.len-1)))) #proportion of male cases in young MSM
  out.p.msm.o.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.msm1+msm.len):(prev.cf.p.diag.msm1+msm.len*2-1)))) #proportion of male cases in old MSM
  out.p.hiv.y.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.p.diag.hiv1:(prev.cf.p.diag.hiv1+hiv.len-1)))) #proportion of young MSM cases with HIV
  out.p.hiv.o.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.hiv1+hiv.len):(prev.cf.p.diag.hiv1+hiv.len*2-1)))) #proportion of old MSM cases with HIV
  out.p.sec.y.m.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.p.diag.sec1:(prev.cf.p.diag.sec1+diag.len-1)))) #proportion of early cases that are secondary
  out.p.sec.o.m.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.sec1+diag.len):(prev.cf.p.diag.sec1+diag.len*2-1)))) #proportion of early cases that are secondary
  out.p.sec.y.f.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.sec1+diag.len*2):(prev.cf.p.diag.sec1+diag.len*3-1)))) #proportion of early cases that are secondary
  out.p.el.y.m.cf <- reshape::melt(as.matrix(subset(pred, select=prev.cf.p.diag.el1:(prev.cf.p.diag.el1+diag.len-1)))) #proportion of early cases that are early latent
  out.p.el.o.m.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.el1+diag.len):(prev.cf.p.diag.el1+diag.len*2-1)))) #proportion of early cases that are early latent
  out.p.el.y.f.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.el1+diag.len*2):(prev.cf.p.diag.el1+diag.len*3-1)))) #proportion of early cases that are early latent
  out.p.early.cf <- reshape::melt(as.matrix(subset(pred, select=(prev.cf.p.diag.early)))) #proportion of early cases that are early latent
  }
  out.inc.y.m.1 <-reshape::melt(100*as.matrix(subset(pred, select=prev.n.inc1:(prev.n.inc1+cal.period-1)))/sum(n.i[y.m.1])) #model incidence by age, sex, and subpopulation
  out.inc.o.m.1 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period):(prev.n.inc1+cal.period*2-1)))/sum(n.i[o.m.1]))
  out.inc.y.m.2 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*2):(prev.n.inc1+cal.period*3-1)))/sum(n.i[y.m.2]))
  out.inc.o.m.2 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*3):(prev.n.inc1+cal.period*4-1)))/sum(n.i[o.m.2]))
  out.inc.y.m.3 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*4):(prev.n.inc1+cal.period*5-1)))/sum(n.i[y.m.3]))
  out.inc.o.m.3 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*5):(prev.n.inc1+cal.period*6-1)))/sum(n.i[o.m.3]))
  out.inc.y.msm <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*6):(prev.n.inc1+cal.period*7-1)))/sum(n.i[y.m.4]))
  out.inc.o.msm <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*7):(prev.n.inc1+cal.period*8-1)))/sum(n.i[o.m.4]))
  out.inc.y.msmhiv <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*8):(prev.n.inc1+cal.period*9-1)))/sum(n.i[y.m.5]))
  out.inc.o.msmhiv <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*9):(prev.n.inc1+cal.period*10-1)))/sum(n.i[o.m.5]))
  out.inc.y.f.1 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*10):(prev.n.inc1+cal.period*11-1)))/sum(n.i[y.f.1]))
  out.inc.o.f.1 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*11):(prev.n.inc1+cal.period*12-1)))/sum(n.i[o.f.1]))
  out.inc.y.f.2 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*12):(prev.n.inc1+cal.period*13-1)))/sum(n.i[y.f.2]))
  out.inc.o.f.2 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*13):(prev.n.inc1+cal.period*14-1)))/sum(n.i[o.f.2]))
  out.inc.y.f.3 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*14):(prev.n.inc1+cal.period*15-1)))/sum(n.i[y.f.3]))
  out.inc.o.f.3 <-reshape::melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*15):(prev.n.inc1+cal.period*16-1)))/sum(n.i[o.f.3]))
  if(exists("showCounterfactual") && showCounterfactual == TRUE){
  out.inc.y.m.1.cf <-reshape::melt(100*as.matrix(subset(pred, select=prev.cf.n.inc1:(prev.cf.n.inc1+cal.period-1)))/sum(n.i[y.m.1])) #model incidence by age, sex, and subpopulation
  out.inc.o.m.1.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period):(prev.cf.n.inc1+cal.period*2-1)))/sum(n.i[o.m.1]))
  out.inc.y.m.2.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*2):(prev.cf.n.inc1+cal.period*3-1)))/sum(n.i[y.m.2]))
  out.inc.o.m.2.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*3):(prev.cf.n.inc1+cal.period*4-1)))/sum(n.i[o.m.2]))
  out.inc.y.m.3.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*4):(prev.cf.n.inc1+cal.period*5-1)))/sum(n.i[y.m.3]))
  out.inc.o.m.3.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*5):(prev.cf.n.inc1+cal.period*6-1)))/sum(n.i[o.m.3]))
  out.inc.y.msm.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*6):(prev.cf.n.inc1+cal.period*7-1)))/sum(n.i[y.m.4]))
  out.inc.o.msm.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*7):(prev.cf.n.inc1+cal.period*8-1)))/sum(n.i[o.m.4]))
  out.inc.y.msmhiv.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*8):(prev.cf.n.inc1+cal.period*9-1)))/sum(n.i[y.m.5]))
  out.inc.o.msmhiv.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*9):(prev.cf.n.inc1+cal.period*10-1)))/sum(n.i[o.m.5]))
  out.inc.y.f.1.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*10):(prev.cf.n.inc1+cal.period*11-1)))/sum(n.i[y.f.1]))
  out.inc.o.f.1.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*11):(prev.cf.n.inc1+cal.period*12-1)))/sum(n.i[o.f.1]))
  out.inc.y.f.2.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*12):(prev.cf.n.inc1+cal.period*13-1)))/sum(n.i[y.f.2]))
  out.inc.o.f.2.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*13):(prev.cf.n.inc1+cal.period*14-1)))/sum(n.i[o.f.2]))
  out.inc.y.f.3.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*14):(prev.cf.n.inc1+cal.period*15-1)))/sum(n.i[y.f.3]))
  out.inc.o.f.3.cf <-reshape::melt(100*as.matrix(subset(pred, select=(prev.cf.n.inc1+cal.period*15):(prev.cf.n.inc1+cal.period*16-1)))/sum(n.i[o.f.3]))
  }
  #browser()
  # For any data separated into compartments, change data frame X2 value to avoid the bug where
  # X2 value of [name]10 gets put before an X2 value of [name][single digit > 1] in the resulting plot.
  # This bug only affects display of older male populations, but in case we change the order of
  # compartments or something, it seems appropriate to do it for all comparents in the affected data sets.
  out.inc.y.m.1$X2 <- sapply(out.inc.y.m.1$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.1$X2 <- sapply(out.inc.o.m.1$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.m.2$X2 <- sapply(out.inc.y.m.2$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.2$X2 <- sapply(out.inc.o.m.2$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.m.3$X2 <- sapply(out.inc.y.m.3$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.3$X2 <- sapply(out.inc.o.m.3$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.msm$X2 <- sapply(out.inc.y.msm$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.msm$X2 <- sapply(out.inc.o.msm$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.msmhiv$X2 <- sapply(out.inc.y.msmhiv$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.msmhiv$X2 <- sapply(out.inc.o.msmhiv$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.1$X2 <- sapply(out.inc.y.f.1$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.1$X2 <- sapply(out.inc.o.f.1$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.2$X2 <- sapply(out.inc.y.f.2$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.2$X2 <- sapply(out.inc.o.f.2$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.3$X2 <- sapply(out.inc.y.f.3$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.3$X2 <- sapply(out.inc.o.f.3$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  out.inc.y.m.1.cf$X2 <- sapply(out.inc.y.m.1.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.1.cf$X2 <- sapply(out.inc.o.m.1.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.m.2.cf$X2 <- sapply(out.inc.y.m.2.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.2.cf$X2 <- sapply(out.inc.o.m.2.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.m.3.cf$X2 <- sapply(out.inc.y.m.3.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.m.3.cf$X2 <- sapply(out.inc.o.m.3.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.msm.cf$X2 <- sapply(out.inc.y.msm.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.msm.cf$X2 <- sapply(out.inc.o.msm.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.msmhiv.cf$X2 <- sapply(out.inc.y.msmhiv.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.msmhiv.cf$X2 <- sapply(out.inc.o.msmhiv.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.1.cf$X2 <- sapply(out.inc.y.f.1.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.1.cf$X2 <- sapply(out.inc.o.f.1.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.2.cf$X2 <- sapply(out.inc.y.f.2.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.2.cf$X2 <- sapply(out.inc.o.f.2.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.y.f.3.cf$X2 <- sapply(out.inc.y.f.3.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.inc.o.f.3.cf$X2 <- sapply(out.inc.o.f.3.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  }
  out.diag.y.m$X2 <- sapply(out.diag.y.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diag.o.m$X2 <- sapply(out.diag.o.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diag.y.f$X2 <- sapply(out.diag.y.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x)) 
  out.diag.o.f$X2 <- sapply(out.diag.o.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x)) 
  out.diaglate.y.m$X2 <- sapply(out.diaglate.y.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.o.m$X2 <- sapply(out.diaglate.o.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.y.f$X2 <- sapply(out.diaglate.y.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.o.f$X2 <- sapply(out.diaglate.o.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x)) 
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  out.diag.y.m.cf$X2 <- sapply(out.diag.y.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diag.o.m.cf$X2 <- sapply(out.diag.o.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diag.y.f.cf$X2 <- sapply(out.diag.y.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x)) 
  out.diag.o.f.cf$X2 <- sapply(out.diag.o.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x)) 
  out.diaglate.y.m.cf$X2 <- sapply(out.diaglate.y.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.o.m.cf$X2 <- sapply(out.diaglate.o.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.y.f.cf$X2 <- sapply(out.diaglate.y.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.diaglate.o.f.cf$X2 <- sapply(out.diaglate.o.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))  
  }
  out.p.msm.y$X2 <- sapply(out.p.msm.y$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.msm.o$X2 <- sapply(out.p.msm.o$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.hiv.y$X2 <- sapply(out.p.hiv.y$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.hiv.o$X2 <- sapply(out.p.hiv.o$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.y.m$X2 <- sapply(out.p.sec.y.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.o.m$X2 <- sapply(out.p.sec.o.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.y.f$X2 <- sapply(out.p.sec.y.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.y.m$X2 <- sapply(out.p.el.y.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.o.m$X2 <- sapply(out.p.el.o.m$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.y.f$X2 <- sapply(out.p.el.y.f$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  out.p.msm.y.cf$X2 <- sapply(out.p.msm.y.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.msm.o.cf$X2 <- sapply(out.p.msm.o.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.hiv.y.cf$X2 <- sapply(out.p.hiv.y.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.hiv.o.cf$X2 <- sapply(out.p.hiv.o.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.y.m.cf$X2 <- sapply(out.p.sec.y.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.o.m.cf$X2 <- sapply(out.p.sec.o.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.sec.y.f.cf$X2 <- sapply(out.p.sec.y.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.y.m.cf$X2 <- sapply(out.p.el.y.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.o.m.cf$X2 <- sapply(out.p.el.o.m.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  out.p.el.y.f.cf$X2 <- sapply(out.p.el.y.f.cf$X2, function(x) gsub("([a-zA-Z.])([0-9])$","\\10\\2",x))
  }
  #prep input data for plotting
  diag.data <- as.data.frame(cbind(year=seq((end.year-(nrow(diag.age.sex.dat))+1),end.year,1), rate=as.numeric(c(diag.age.sex.rate.m, diag.age.sex.rate.f.y, diag.age.sex.rate.f.o)), cat=rep(c("y.m", "o.m", "y.f", "o.f"),each=nrow(diag.age.sex.dat))))
  age.data <- as.data.frame(age.dist.dat[1:4,])
  
  ### plot input data and model outputs ###
  
  #get max values of reported cases to set y-axis
  max.y.m <- max(out.diag.y.m$value, na.rm=TRUE)
  max.o.m <- max(out.diag.o.m$value, na.rm=TRUE)
  max.y.f <-max(out.diag.y.f$value, na.rm=TRUE)
  max.o.f <-max(out.diag.o.f$value, na.rm=TRUE)
  max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f,diag.data$diag.age.sex.rate) +5
	max.diag.y <- max(max.y.m, max.y.f) + 5
	max.diag.o <- max(max.o.m, max.o.f) + 5
  
  ### plot reported cases by age and sex ###
  plot.diag.y.m <- ggplot(data=out.diag.y.m)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_m), color="red", shape=15, size=2) +
    labs(title="All M 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    expand_limits(y=c(0, max.diag.y))+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.o.m <- ggplot(data=out.diag.o.m)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_m), color="red", shape=15, size=2) +
    labs(title="All M 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    expand_limits(y=c(0, max.diag.y))+
    #expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.y.f <- ggplot(data=out.diag.y.f)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_f), color="red", shape=15, size=2) +
    labs(title="All F 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    #expand_limits(y=0)+
    expand_limits(y=c(0, max.diag.o))+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.o.f <- ggplot(data=out.diag.o.f)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_f), color="red", shape=15, size=2) +
    labs(title="All F 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    expand_limits(y=c(0, max.diag.o))+
    ##expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  plot.diag.y.m.cf <- ggplot(data=out.diag.y.m.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_m), color="red", shape=15, size=2) +
    labs(title="All M 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    #expand_limits(y=0)+
    expand_limits(y=c(0, max.diag))+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.o.m.cf <- ggplot(data=out.diag.o.m.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_m), color="red", shape=15, size=2) +
    labs(title="All M 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    expand_limits(y=c(0, max.diag))+
    #expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.y.f.cf <- ggplot(data=out.diag.y.f.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_f), color="red", shape=15, size=2) +
    labs(title="All F 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    #expand_limits(y=0)+
    expand_limits(y=c(0, max.diag))+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  
  plot.diag.o.f.cf <- ggplot(data=out.diag.o.f.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_f), color="red", shape=15, size=2) +
    labs(title="All F 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
    # coord_cartesian(ylim=c(0,max.diag))+
    expand_limits(y=c(0, max.diag))+
    ##expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))
  }
  ### additional plots ###
  
  plot.rr.diag <- ggplot(data=out.rr.diag)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_pointrange(data=diag.rr, aes(x=1:8, y=mean, ymin=min, ymax=max), color="red", shape=15, size=0.75, lty=3) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Relative diagnosis rates in subpopulations", x="Population", y="Reported case relative risk") +
    scale_x_discrete(labels=c("M black y","M black o", "M Hispanic y", "M Hispanic o", "F black y", "F black o", "F Hispanic y", "F Hispanic o"))
  
  plot.subpop <- ggplot(data=out.subpop)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    #geom_pointrange(data=s.dist.sd,aes(x=1:nrow(s.dist.sd), y=dat.s.dist,ymin=lcl,ymax=ucl), color="red", shape=15, size=0.5) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="Subpopulation assortative mixing", x="Population", y="Proportion of partnerships \nwith same subpopulation") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=c("M black","M other","M Hispanic","F black","F other", "F Hispanic", "MSM HIV-", "MSM HIV+"))
  
  plot.p.msm.y <- ggplot(data=out.p.msm.y)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_point(data=as.data.frame(p.msm.cases[1:msm.len]),aes(x=1:msm.len, y=p.msm.cases[1:msm.len]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="MSM 20-44 y", x="Year", y="Proportion of early syphilis\ncases in MSM") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-msm.len+1,end.year,1))
  
  plot.p.msm.o <- ggplot(data=out.p.msm.o)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_point(data=as.data.frame(p.msm.cases[(msm.len+1):(2*msm.len)]),aes(x=1:msm.len, y=p.msm.cases[(msm.len+1):(2*msm.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="MSM 45-64 y", x="Year", y="Proportion of early syphilis\ncases in MSM") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-msm.len+1,end.year,1))
  
  plot.p.hiv.y <- ggplot(data=out.p.hiv.y)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_point(data=as.data.frame(p.hiv.cases[1:hiv.len]),aes(x=1:hiv.len, y=p.hiv.cases[1:hiv.len]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="MSM 20-44 y", x="Year", y="Proportion of MSM cases\nwith HIV coinfection") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-hiv.len+1,end.year,1))
  
  plot.p.hiv.o <- ggplot(data=out.p.hiv.o)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_point(data=as.data.frame(p.hiv.cases[(hiv.len+1):(2*hiv.len)]),aes(x=1:hiv.len, y=p.hiv.cases[(hiv.len+1):(2*hiv.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="MSM 45-64 y", x="Year", y="Proportion of MSM cases\nwith HIV coinfection") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-hiv.len+1,end.year,1))
  
  plot.p.sec.y.m <- ggplot(data=out.p.sec.y.m)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.sec[1:diag.len]),aes(x=1:diag.len, y=p.sec[1:diag.len]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="M 20-44 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.p.sec.o.m <- ggplot(data=out.p.sec.o.m)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.sec[(diag.len+1):(2*diag.len)]),aes(x=1:diag.len, y=p.sec[(diag.len+1):(2*diag.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="M 45-64 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.p.sec.y.f <- ggplot(data=out.p.sec.y.f)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.sec[(2*diag.len+1):(3*diag.len)]),aes(x=1:diag.len, y=p.sec[(2*diag.len+1):(3*diag.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="F 20-44 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.p.el.y.m <- ggplot(data=out.p.el.y.m)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.el[1:diag.len]),aes(x=1:diag.len, y=p.el[1:diag.len]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="M 20-44 y", x="Year", y="Proportion of reported early cases\nthat are early latent") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.p.el.o.m <- ggplot(data=out.p.el.o.m)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.el[(diag.len+1):(2*diag.len)]),aes(x=1:diag.len, y=p.el[(diag.len+1):(2*diag.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="M 45-64 y", x="Year", y="Proportion of diagnosed early cases\nthat are early latent") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.p.el.y.f <- ggplot(data=out.p.el.y.f)+
    geom_boxplot(aes(x=X2, y=value), fill="purple1") +
    geom_point(data=as.data.frame(p.el[(2*diag.len+1):(3*diag.len)]),aes(x=1:diag.len, y=p.el[(2*diag.len+1):(3*diag.len)]), color="red", shape=15, size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8) ) + 
    labs(title="F 20-44 y", x="Year", y="Proportion of diagnosed early cases\nthat are early latent") +
    coord_cartesian(ylim=c(0,1))+
    scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))
  
  plot.age <- ggplot(data=out.age)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_pointrange(data=age.data[1,],aes(x=1, y=p.same.age, ymin=lcl, ymax=ucl), color="red", shape=15, size=0.5) +
    geom_pointrange(data=age.data[2,],aes(x=3, y=p.same.age, ymin=lcl, ymax=ucl), color="red", shape=15, size=0.5) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8)) + 
    labs(title="Age assortative mixing", x="Population", y="Proportion of partnerships \nwith same age group") +
    coord_cartesian(ylim=c(0, 1))+
    scale_x_discrete(labels=c("M young","M old","F young", "F old"))
  
  plot.p.early <- ggplot(data=out.p.early)+
    geom_boxplot(aes(x=X2, y=value, fill=X2)) +
    geom_pointrange(data=p.early.dat,aes(x=1, y=mean, ymin=min, ymax=max), color="red", shape=15, size=0.5) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1),title=element_text(size=8)) + 
    labs(title="Proportion early syphilis cases",x="", y="Proportion of all reported syphilis cases \nthat are primary, secondary, or early latent") +
    coord_cartesian(ylim=c(0, 1)) +
    scale_x_discrete(labels="")
  
  
  
  ### plot incident cases by subpop, age, and sex ###
  #get max values of incidence for setting y-axis 
  max.inc.y.m <- max(out.inc.y.m.1$value, out.inc.y.m.2$value, out.inc.y.m.3$value, na.rm=TRUE)
  max.inc.o.m <- max(out.inc.o.m.1$value, out.inc.o.m.2$value, out.inc.o.m.3$value, na.rm=TRUE)
  max.inc.y.f <-max(out.inc.y.f.1$value, out.inc.y.f.2$value, out.inc.y.f.3$value, na.rm=TRUE)
  max.inc.o.f <-max(out.inc.o.f.1$value, out.inc.o.f.2$value, out.inc.o.f.3$value, na.rm=TRUE) 
  max.inc <- ceiling(max(max.inc.y.m, max.inc.o.m, max.inc.y.f, max.inc.o.f))
  max.inc.msm <-ceiling(max(out.inc.y.msm$value, out.inc.o.msm$value, na.rm=TRUE)+1)
  ### plot incident cases by subpop, age, and sex: 2000-2015 ###
  
  plot.inc.y.m.1 <- ggplot(data=out.inc.y.m.1)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.1 <- ggplot(data=out.inc.o.m.1)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    #scale_x_discrete(labels=seq((end.year-(cal.period*3)+1),end.year,1))
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.m.2 <- ggplot(data=out.inc.y.m.2)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.2<- ggplot(data=out.inc.o.m.2)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.m.3 <- ggplot(data=out.inc.y.m.3)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.3<- ggplot(data=out.inc.o.m.3)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.1 <- ggplot(data=out.inc.y.f.1)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.1 <- ggplot(data=out.inc.o.f.1)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.2 <- ggplot(data=out.inc.y.f.2)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.2 <- ggplot(data=out.inc.o.f.2)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.3 <- ggplot(data=out.inc.y.f.3)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.3 <- ggplot(data=out.inc.o.f.3)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.msm <- ggplot(data=out.inc.y.msm)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV- MSM 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.msm <- ggplot(data=out.inc.o.msm)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV- MSM 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.msmhiv <- ggplot(data=out.inc.y.msmhiv)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV+ MSM 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.msmhiv <- ggplot(data=out.inc.o.msmhiv)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV+ MSM 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  # for the no-contact-tracing counterfactual
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  plot.inc.y.m.1.cf <- ggplot(data=out.inc.y.m.1.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.1.cf <- ggplot(data=out.inc.o.m.1.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    #scale_x_discrete(labels=seq((end.year-(cal.period*3)+1),end.year,1))
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.m.2.cf <- ggplot(data=out.inc.y.m.2.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.2.cf<- ggplot(data=out.inc.o.m.2.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.m.3.cf <- ggplot(data=out.inc.y.m.3.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic M 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.m.3.cf<- ggplot(data=out.inc.o.m.3.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic M 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.1.cf <- ggplot(data=out.inc.y.f.1.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.1.cf <- ggplot(data=out.inc.o.f.1.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Black F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.2.cf <- ggplot(data=out.inc.y.f.2.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.2.cf <- ggplot(data=out.inc.o.f.2.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Other F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.f.3.cf <- ggplot(data=out.inc.y.f.3.cf)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic F 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.f.3.cf <- ggplot(data=out.inc.o.f.3)+
    geom_line(aes(x=X2, y=value, group=X1),color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="Hispanic F 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.msm.cf <- ggplot(data=out.inc.y.msm.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV- MSM 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.msm.cf <- ggplot(data=out.inc.o.msm.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV- MSM 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.y.msmhiv.cf <- ggplot(data=out.inc.y.msmhiv.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV+ MSM 20-44 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  
  plot.inc.o.msmhiv.cf <- ggplot(data=out.inc.o.msmhiv.cf)+
    geom_line(aes(x=X2, y=value, group=X1), color="grey") +
    stat_summary(aes(x=X2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
    theme(legend.position="none", axis.text.x=element_text(angle=90,hjust=1), title=element_text(size=8)) + 
    labs(title="HIV+ MSM 45-64 y", x="Year", y="Incidence (%)") +
    #coord_cartesian(ylim=c(0,max.inc.msm))+
    expand_limits(y=0)+
    scale_x_discrete(labels=seq((end.year-cal.period+1),end.year,1))
  }
  
  #######################################  
  ####   Plot priors and posteriors   ###
  #######################################      
  
  #epsilon.1
  epsilon.1.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.epsilon.1"])))
  epsilon.1.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.1"],prior.param2["epsilon.1"])))
  plot.epsilon.1 <- ggplot() +
    geom_histogram(data=epsilon.1.post, aes(x=X1, y=..density..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Risk group mixing:\n black")
  
  #epsilon.2
  epsilon.2.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.epsilon.2"])))
  epsilon.2.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.2"],prior.param2["epsilon.2"])))
  plot.epsilon.2 <- ggplot() +
    geom_histogram(data=epsilon.2.post, aes(x=X1, y=..density..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Risk group mixing:\n other")
  
  #epsilon.3
  epsilon.3.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.epsilon.3"])))
  epsilon.3.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.3"],prior.param2["epsilon.3"])))
  plot.epsilon.3 <- ggplot() +
    geom_histogram(data=epsilon.3.post, aes(x=X1, y=..density..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.3.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Risk group mixing:\n Hispanic")
  
  #epsilon.4
  epsilon.4.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.epsilon.4"])))
  epsilon.4.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.4"],prior.param2["epsilon.4"])))
  plot.epsilon.4 <- ggplot() +
    geom_histogram(data=epsilon.4.post, aes(x=X1, y=..density..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.4.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Risk group mixing:\n HIV- MSM")
  
  #epsilon.5
  epsilon.5.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.epsilon.5"])))
  epsilon.5.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["epsilon.5"],prior.param2["epsilon.5"])))
  plot.epsilon.5 <- ggplot() +
    geom_histogram(data=epsilon.5.post, aes(x=X1, y=..density..), fill="slateblue",size=0.1, colour="slateblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=epsilon.5.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Risk group mixing:\n HIV+ MSM")
  
  #pi.m
  pi.m.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.pi.m"])))
  pi.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.m"],prior.param2["pi.m"])))
  plot.pi.m <- ggplot() +
    geom_histogram(data=pi.m.post, aes(x=X1, y=..density..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Age assortativity:\n young M/old F")
  
  #pi.f
  pi.f.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.pi.f"])))
  pi.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.f"],prior.param2["pi.f"])))
  plot.pi.f <- ggplot() +
    geom_histogram(data=pi.f.post, aes(x=X1, y=..density..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Age assortativity:\n young F/old M")
  
  #pi.msm
  pi.msm.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.pi.msm"])))
  pi.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["pi.msm"],prior.param2["pi.msm"])))
  plot.pi.msm <- ggplot() +
    geom_histogram(data=pi.msm.post, aes(x=X1, y=..density..), fill="dodgerblue",size=0.1, colour="dodgerblue", alpha=0.6, binwidth=0.025)+
    geom_area(data=pi.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x="Age assortativity:\n MSM")
  
  #theta.1
  theta.1.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.1"])))
  theta.1.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.1"],prior.param2["theta.1"])))
  plot.theta.1 <- ggplot() +
    geom_histogram(data=theta.1.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n black M")
  
  #theta.2
  theta.2.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.2"])))
  theta.2.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.2"],prior.param2["theta.2"])))
  plot.theta.2 <- ggplot() +
    geom_histogram(data=theta.2.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n other M")
  
  #theta.3
  theta.3.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.3"])))
  theta.3.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.3"],prior.param2["theta.3"])))
  plot.theta.3 <- ggplot() +
    geom_histogram(data=theta.3.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.3.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n Hispanic M")
  
  #theta.4
  theta.4.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.4"])))
  theta.4.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.4"],prior.param2["theta.4"])))
  plot.theta.4 <- ggplot() +
    geom_histogram(data=theta.4.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.4.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n HIV- MSM")
  
  #theta.5
  theta.5.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.5"])))
  theta.5.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.5"],prior.param2["theta.5"])))
  plot.theta.5 <- ggplot() +
    geom_histogram(data=theta.5.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.5.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n HIV+ MSM")
  
  #theta.6
  theta.6.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.6"])))
  theta.6.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.6"],prior.param2["theta.6"])))
  plot.theta.6 <- ggplot() +
    geom_histogram(data=theta.6.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.6.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n black F")
  
  #theta.7
  theta.7.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.7"])))
  theta.7.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.7"],prior.param2["theta.7"])))
  plot.theta.7 <- ggplot() +
    geom_histogram(data=theta.7.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.7.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n other F")
  
  #theta.8
  theta.8.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.8"])))
  theta.8.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.8"],prior.param2["theta.8"])))
  plot.theta.8 <- ggplot() +
    geom_histogram(data=theta.8.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.8.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.5,1))+
    labs(x="Subpopulation mixing:\n Hispanic F")
  
  #theta.hiv
  theta.hiv.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.theta.hiv"])))
  theta.hiv.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["theta.hiv"],prior.param2["theta.hiv"])))
  plot.theta.hiv <- ggplot() +
    geom_histogram(data=theta.hiv.post, aes(x=X1, y=..density..), fill="cornflowerblue",size=0.1, colour="cornflowerblue", alpha=0.6, binwidth=0.01)+
    geom_area(data=theta.hiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Mixing with MSM of same HIV status")))
  
  #c.min.m1
  c.min.m1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.m.1"])))
  c.min.m1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.m1"],prior.param2["c.min.m1"])))
  plot.c.min.m1 <- ggplot() +
    geom_histogram(data=c.min.m1.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.m1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " M 20-44y")))
  
  #c.min.m2
  c.min.m2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.m.2"])))
  c.min.m2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.m2"],prior.param2["c.min.m2"])))
  plot.c.min.m2 <- ggplot() +
    geom_histogram(data=c.min.m2.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.m2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " M 45-64 y")))
  
  #c.min.f1
  c.min.f1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.f.1"])))
  c.min.f1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.f1"],prior.param2["c.min.f1"])))
  plot.c.min.f1 <- ggplot() +
    geom_histogram(data=c.min.f1.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.f1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " F 20-44 y")))
  
  #c.min.f2
  c.min.f2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.f.2"])))
  c.min.f2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.f2"],prior.param2["c.min.f2"])))
  plot.c.min.f2 <- ggplot() +
    geom_histogram(data=c.min.f2.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.f2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " F 45-64y")))
  
  #c.min.msm1
  c.min.msm1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.msm.1"])))
  c.min.msm1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.msm1"],prior.param2["c.min.msm1"])))
  plot.c.min.msm1 <- ggplot() +
    geom_histogram(data=c.min.msm1.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.msm1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " MSM 20-44 y")))
  
  #c.min.msm2
  c.min.msm2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.c.min.msm.2"])))
  c.min.msm2.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["c.min.msm2"],prior.param2["c.min.msm2"])))
  plot.c.min.msm2 <- ggplot() +
    geom_histogram(data=c.min.msm2.post, aes(x=X1, y=..density..), fill="lightgoldenrod1",size=0.1, colour="lightgoldenrod3", alpha=0.6, binwidth=0.05)+
    geom_area(data=c.min.msm2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("c"[min], " MSM 45-64 y")))
  
  #dur.incub
  dur.incub.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.incub"])))
  dur.incub.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dnorm(seq(0,50,0.1),prior.param1["dur.incub"],prior.param2["dur.incub"])))
  plot.dur.incub <- ggplot() +
    geom_histogram(data=dur.incub.post, aes(x=X1, y=..density..), fill="turquoise",size=0.5, colour="turquoise", alpha=0.6, binwidth=1)+
    geom_area(data=dur.incub.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,40))+
    labs(x=expression(paste("Latent period (d)")))
  
  #dur.prim
  dur.prim.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.prim"])))
  dur.prim.prior <- as.data.frame(cbind(x=seq(0,75,0.1),y=dnorm(seq(0,75,0.1),prior.param1["dur.prim"],prior.param2["dur.prim"])))
  plot.dur.prim <- ggplot() +
    geom_histogram(data=dur.prim.post, aes(x=X1, y=..density..), fill="turquoise",size=0.5, colour="turquoise", alpha=0.6, binwidth=2)+
    geom_area(data=dur.prim.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,75))+
    labs(x=expression(paste("Duration, primary (d)")))
  
  #dur.sec
  dur.sec.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.sec"])))
  dur.sec.prior <- as.data.frame(cbind(x=seq(0,150,0.1),y=dunif(seq(0,150,0.1),prior.param1["dur.sec"],prior.param2["dur.sec"])))
  plot.dur.sec <- ggplot() +
    geom_histogram(data=dur.sec.post, aes(x=X1, y=..density..), fill="turquoise",size=0.5, colour="turquoise", alpha=0.6, binwidth=2)+
    geom_area(data=dur.sec.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,150))+
    labs(x=expression(paste("Duration, secondary (d)")))
  
  #dur.imm.inf
  dur.imm.inf.post<- data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.imm.inf"])))
  dur.imm.inf.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["dur.imm.inf"],prior.param2["dur.imm.inf"])))
  plot.dur.imm.inf <<- ggplot() +
    geom_histogram(data=dur.imm.inf.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=1)+
    geom_area(data=dur.imm.inf.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("Dur immune, treated P&S (d)")))
  
  #dur.imm.early
  dur.imm.early.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.imm.early"])))
  dur.imm.early.prior <- as.data.frame(cbind(x=seq(0,20,0.1),y=dunif(seq(0,20,0.1),prior.param1["dur.imm.early"],prior.param2["dur.imm.early"])))
  plot.dur.imm.early <<- ggplot() +
    geom_histogram(data=dur.imm.early.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.1)+
    geom_area(data=dur.imm.early.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("Dur immune multiplier, treated early latent")))
  
  #dur.immune
  dur.immune.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.dur.immune"])))
  dur.immune.prior <- as.data.frame(cbind(x=seq(0,3000,1),y=dnorm(seq(0,3000,1),prior.param1["dur.immune"],prior.param2["dur.immune"])))
  plot.dur.immune <<- ggplot() +
    geom_histogram(data=dur.immune.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=20)+
    geom_area(data=dur.immune.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=5))) + 
    coord_cartesian(xlim=c(0,2000))+
    labs(x=expression(paste("Dur immune, treated late latent (d)")))
  
  #b.m
  b.m.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.b.m"])))
  b.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.m"],prior.param2["b.m"])))
  plot.b.m <- ggplot() +
    geom_histogram(data=b.m.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    geom_area(data=b.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Transmission probability:\n F to M")))
  
  #b.f
  b.f.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.b.f"])))
  b.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.f"],prior.param2["b.f"])))
  plot.b.f <- ggplot() +
    geom_histogram(data=b.f.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    geom_area(data=b.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    theme(axis.title.x=element_text(size=8))+
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Transmission probability:\nM to F")))
  
  #b.msm
  b.msm.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.b.msm"])))
  b.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["b.msm"],prior.param2["b.msm"])))
  plot.b.msm <- ggplot() +
    geom_histogram(data=b.msm.post, aes(x=X1, y=..density..), fill="dodgerblue3",size=0.5, colour="dodgerblue3", alpha=0.6, binwidth=0.025)+
    geom_area(data=b.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Transmission probability:\nM to M")))
  
  #rr.screen.m1
  rr.screen.m1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.m1"])))
  rr.screen.m1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.m1"],prior.param2["rr.screen.m1"])))
  plot.rr.screen.m1 <- ggplot() +
    geom_histogram(data=rr.screen.m1.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.m1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nBlack M")))
  
  #rr.screen.m3
  rr.screen.m3.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.m3"])))
  rr.screen.m3.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.m3"],prior.param2["rr.screen.m3"])))
  plot.rr.screen.m3 <- ggplot() +
    geom_histogram(data=rr.screen.m3.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.m3.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nHispanic M")))
  
  #rr.screen.o.m
  rr.screen.o.m.post <- data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.o.m"])))
  rr.screen.o.m.prior <- as.data.frame(cbind(x=seq(0,20,0.01),y=dgamma(seq(0,20,0.01),prior.param1["rr.screen.o.m"],prior.param2["rr.screen.o.m"])))
  plot.rr.screen.o.m <- ggplot() +
    geom_histogram(data=rr.screen.o.m.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.o.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,20))+
    labs(x=expression(paste("RR screen:\nM 45-64y")))
  
  #rr.screen.msmhiv
  rr.screen.msmhiv.post<- data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.msmhiv"])))
  rr.screen.msmhiv.prior <- as.data.frame(cbind(x=seq(0,20,0.01),y=dgamma(seq(0,20,0.01),prior.param1["rr.screen.msmhiv"],prior.param2["rr.screen.msmhiv"])))
  plot.rr.screen.msmhiv <- ggplot() +
    geom_histogram(data=rr.screen.msmhiv.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.msmhiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,20))+
    labs(x=expression(paste("RR screen:\nHIV+ MSM")))
  
  #rr.screen.o.msm
  rr.screen.o.msm.post <- data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.o.msm"])))
  rr.screen.o.msm.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.o.msm"],prior.param2["rr.screen.o.msm"])))
  plot.rr.screen.o.msm <- ggplot() +
    geom_histogram(data=rr.screen.o.msm.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.o.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nMSM 45-64y")))
  
  ##rr.screen.f1
  rr.screen.f1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.f1"])))
  rr.screen.f1.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.f1"],prior.param2["rr.screen.f1"])))
  plot.rr.screen.f1 <- ggplot() +
    geom_histogram(data=rr.screen.f1.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.05)+
    geom_area(data=rr.screen.f1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +  
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nblack F")))
  
  ##rr.screen.f3
  rr.screen.f3.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.f3"])))
  rr.screen.f3.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.f3"],prior.param2["rr.screen.f3"])))
  plot.rr.screen.f3 <- ggplot() +
    geom_histogram(data=rr.screen.f3.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.f3.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nHispanic F")))
  
  ##rr.screen.o.f
  rr.screen.o.f.post <- data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.o.f"])))
  rr.screen.o.f.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.o.f"],prior.param2["rr.screen.o.f"])))
  plot.rr.screen.o.f <- ggplot() +
    geom_histogram(data=rr.screen.o.f.post, aes(x=X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.o.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\n F 45-64y")))
  
  ##rr.screen.ac
  rr.screen.ac.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rr.screen.ac"])))
  rr.screen.ac.prior <- as.data.frame(cbind(x=seq(0,5,0.01),y=dgamma(seq(0,5,0.01),prior.param1["rr.screen.ac"],prior.param2["rr.screen.ac"])))
  plot.rr.screen.ac <- ggplot() +
    geom_histogram(data=rr.screen.ac.post, aes(x=1+X1, y=..density..), fill="deepskyblue",size=0.5, colour="deepskyblue", alpha=0.6, binwidth=0.1)+
    geom_area(data=rr.screen.ac.prior,aes(x=1+x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR screen:\nhigh sexual activity group")))
  
  #p.trt.prim.m
  trt.prim.m.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.prim.m"])))
  trt.prim.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.prim.m"],prior.param2["p.trt.prim.m"])))
  plot.trt.prim.m <- ggplot() +
    geom_histogram(data=trt.prim.m.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.prim.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate primary (/y):\nM")))
  
  #p.trt.prim.f
  trt.prim.f.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.prim.f"])))
  trt.prim.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.prim.f"],prior.param2["p.trt.prim.f"])))
  plot.trt.prim.f <- ggplot() +
    geom_histogram(data=trt.prim.f.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.prim.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate primary (/y):\nF")))
  
  #p.trt.prim.msm
  trt.prim.msm.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.prim.msm"])))
  trt.prim.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.prim.msm"],prior.param2["p.trt.prim.msm"])))
  plot.trt.prim.msm <- ggplot() +
    geom_histogram(data=trt.prim.msm.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.prim.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate primary (/y):\nHIV- MSM ")))
  
  #p.trt.prim.msmhiv
  trt.prim.msmhiv.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.prim.msmhiv"])))
  trt.prim.msmhiv.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.prim.msmhiv"],prior.param2["p.trt.prim.msmhiv"])))
  plot.trt.prim.msmhiv <- ggplot() +
    geom_histogram(data=trt.prim.msmhiv.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.prim.msmhiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate primary (/y):\nHIV+ MSM")))
  
  #p.trt.sec.m
  trt.sec.m.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.sec.m"])))
  trt.sec.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.sec.m"],prior.param2["p.trt.sec.m"])))
  plot.trt.sec.m <- ggplot() +
    geom_histogram(data=trt.sec.m.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.sec.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate secondary (/y):\nM")))
  
  #p.trt.sec.f
  trt.sec.f.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.sec.f"])))
  trt.sec.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.sec.f"],prior.param2["p.trt.sec.f"])))
  plot.trt.sec.f <- ggplot() +
    geom_histogram(data=trt.sec.f.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.sec.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate secondary (/y):F")))
  
  #p.trt.sec.msm
  trt.sec.msm.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.sec.msm"])))
  trt.sec.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.sec.msm"],prior.param2["p.trt.sec.msm"])))
  plot.trt.sec.msm <- ggplot() +
    geom_histogram(data=trt.sec.msm.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.sec.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate secondary (/y):\nHIV- MSM"))) 
  
  #p.trt.sec.msmhiv
  trt.sec.msmhiv.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.sec.msmhiv"])))
  trt.sec.msmhiv.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.sec.msmhiv"],prior.param2["p.trt.sec.msmhiv"])))
  plot.trt.sec.msmhiv <- ggplot() +
    geom_histogram(data=trt.sec.msmhiv.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.sec.msmhiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate secondary (/y):\nHIV+ MSM")))
  
  #p.trt.lat.m
  trt.lat.m.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.lat.m"])))
  trt.lat.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.lat.m"],prior.param2["p.trt.lat.m"])))
  plot.trt.lat.m <- ggplot() +
    geom_histogram(data=trt.lat.m.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.lat.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate early latent (/y):\nM")))
  
  #p.trt.lat.f
  trt.lat.f.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.lat.f"])))
  trt.lat.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.lat.f"],prior.param2["p.trt.lat.f"])))
  plot.trt.lat.f <- ggplot() +
    geom_histogram(data=trt.lat.f.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.lat.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate early latent (/y):\nF")))
  
  #p.trt.lat.msm
  trt.lat.msm.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.lat.msm"])))
  trt.lat.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.lat.msm"],prior.param2["p.trt.lat.msm"])))
  plot.trt.lat.msm <- ggplot() +
    geom_histogram(data=trt.lat.msm.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.lat.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate early latent (/y):\nHIV- MSM")))
  
  #p.trt.lat.msmhiv
  trt.lat.msmhiv.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.lat.msmhiv"])))
  trt.lat.msmhiv.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.lat.msmhiv"],prior.param2["p.trt.lat.msmhiv"])))
  plot.trt.lat.msmhiv <- ggplot() +
    geom_histogram(data=trt.lat.msmhiv.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.lat.msmhiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate early latent (/y):\nHIV+ MSM"))) 
  
  #p.trt.late.m
  trt.late.m.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.late.m"])))
  trt.late.m.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.late.m"],prior.param2["p.trt.late.m"])))
  plot.trt.late.m <- ggplot() +
    geom_histogram(data=trt.late.m.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.late.m.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1)) +
    labs(x=expression(paste("Trt rate late latent (/y):\nM"))) 
  
  #p.trt.late.f
  trt.late.f.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.late.f"])))
  trt.late.f.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.late.f"],prior.param2["p.trt.late.f"])))
  plot.trt.late.f <- ggplot() +
    geom_histogram(data=trt.late.f.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.late.f.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1)) +
    labs(x=expression(paste("Trt rate late latent (/y):\nF"))) 
  
  #p.trt.late.msm
  trt.late.msm.post<- data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.late.msm"])))
  trt.late.msm.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.late.msm"],prior.param2["p.trt.late.msm"])))
  plot.trt.late.msm <- ggplot() +
    geom_histogram(data=trt.late.msm.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.late.msm.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1))+
    labs(x=expression(paste("Trt rate late latent (/y):\nHIV- MSM")))
  
  #p.trt.late.msmhiv
  trt.late.msmhiv.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.p.trt.late.msmhiv"])))
  trt.late.msmhiv.prior <- as.data.frame(cbind(x=seq(0,1,0.01),y=dbeta(seq(0,1,0.01),prior.param1["p.trt.late.msmhiv"],prior.param2["p.trt.late.msmhiv"])))
  plot.trt.late.msmhiv <- ggplot() +
    geom_histogram(data=trt.late.msmhiv.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.01)+
    geom_area(data=trt.late.msmhiv.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,1)) +
    labs(x=expression(paste("Trt rate late latent (/y):\nHIV+ MSM")))
  
  #rr.rep.symp.f
  risk.rep.symp.m.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.rr.rep.symp.m"])))
  risk.rep.symp.m.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["rr.rep.symp.m"],prior.param2["rr.rep.symp.m"])))
  plot.rr.rep.symp.m <- ggplot() +
    geom_histogram(data=risk.rep.symp.m.post, aes(x=X1, y=..density..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.m.prior,aes(x=x, y=y),fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.0,1.1))+
    labs(x=expression(paste("RR report\n M seeking medical treatment")))
  
  #rr.rep.symp.f
  risk.rep.symp.f.post<-data.frame(X1=as.numeric(ilogit(trace.burn.thin[,"logit.rr.rep.symp.f"])))
  risk.rep.symp.f.prior <- as.data.frame(cbind(x=seq(0,1,0.001),y=dbeta(seq(0,1,0.001),prior.param1["rr.rep.symp.f"],prior.param2["rr.rep.symp.f"])))
  plot.rr.rep.symp.f <- ggplot() +
    geom_histogram(data=risk.rep.symp.f.post, aes(x=X1, y=..density..), fill="thistle3",size=0.5, colour="thistle3", alpha=0.6, binwidth=0.05)+
    geom_area(data=risk.rep.symp.f.prior,aes(x=x, y=y),fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8),axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0.0,1.1))+
    labs(x=expression(paste("RR report\nF seeking medical treatment")))
  
  #rp.1.1.1.1
  rp.1.1.1.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.1.1.1"])))
  rp.1.1.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.1.1.1.1"],prior.param2["rp.1.1.1.1"])))
  plot.rp.1.1.1.1 <- ggplot() +
    geom_histogram(data=rp.1.1.1.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.1.1.1.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nBlack M low AC, 20-44 y")))
  
  #rp.1.1.2.1
  rp.1.1.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.1.2.1"])))
  rp.1.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.1.1.2.1"],prior.param2["rp.1.1.2.1"])))
  plot.rp.1.1.2.1 <- ggplot() +
    geom_histogram(data=rp.1.1.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.1.1.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nBlack M high AC, 20-44 y")))
  
  #rp.1.2.1.1
  rp.1.2.1.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.2.1.1"])))
  rp.1.2.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.1.2.1.1"],prior.param2["rp.1.2.1.1"])))
  plot.rp.1.2.1.1 <- ggplot() +
    geom_histogram(data=rp.1.2.1.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.1.2.1.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nBlack F low AC, 20-44 y")))
  
  #rp.1.2.2.1
  rp.1.2.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.2.2.1"])))
  rp.1.2.2.1.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.1.2.2.1"],prior.param2["rp.1.2.2.1"])))
  plot.rp.1.2.2.1 <- ggplot() +
    geom_histogram(data=rp.1.2.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.1.2.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nBlack F high AC, 20-44 y")))
  
  #rp.1.1.1.2
  rp.1.1.1.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.1.1.2"])))
  rp.1.1.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.1.1.1.2"],prior.param2["rp.1.1.1.2"])))
  plot.rp.1.1.1.2 <- ggplot() +
    geom_histogram(data=rp.1.1.1.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.1.1.1.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR partner change:\nBlack M low AC, 45-64 y")))
  
  #rp.1.1.2.2
  rp.1.1.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.1.2.2"])))
  rp.1.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.1.1.2.2"],prior.param2["rp.1.1.2.2"])))
  plot.rp.1.1.2.2 <- ggplot() +
    geom_histogram(data=rp.1.1.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.1.1.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nBlack M high AC, 45-64 y")))
  
  #rp.1.2.1.2
  rp.1.2.1.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.2.1.2"])))
  rp.1.2.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.1.2.1.2"],prior.param2["rp.1.2.1.2"])))
  plot.rp.1.2.1.2 <- ggplot() +
    geom_histogram(data=rp.1.2.1.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.1.2.1.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR partner change:\nBlack F low AC, 45-64 y")))
  
  #rp.1.2.2.2
  rp.1.2.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.2.2.2"])))
  rp.1.2.2.2.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.1.2.2.2"],prior.param2["rp.1.2.2.2"])))
  plot.rp.1.2.2.2 <- ggplot() +
    geom_histogram(data=rp.1.2.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.1.2.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nBlack F high AC, 45-64 y")))
  
  #rp.2.1.1.1
  rp.2.1.1.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.1.1.1"])))
  rp.2.1.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.2.1.1.1"],prior.param2["rp.2.1.1.1"])))
  plot.rp.2.1.1.1 <- ggplot() +
    geom_histogram(data=rp.2.1.1.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.1.1.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nOther M low AC, 20-44 y")))
  
  #rp.2.1.2.1
  rp.2.1.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.1.2.1"])))
  rp.2.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.2.1.2.1"],prior.param2["rp.2.1.2.1"])))
  plot.rp.2.1.2.1 <- ggplot() +
    geom_histogram(data=rp.2.1.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.2.1.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nOther M high AC, 20-44 y")))
  
  #rp.2.2.1.1
  rp.2.2.1.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.2.1.1"])))
  rp.2.2.1.1.prior <- as.data.frame(cbind(x=seq(0,10,0.1),y=dgamma(seq(0,10,0.1),prior.param1["rp.2.2.1.1"],prior.param2["rp.2.2.1.1"])))
  plot.rp.2.2.1.1 <- ggplot() +
    geom_histogram(data=rp.1.2.1.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.5)+
    geom_area(data=rp.2.2.1.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nOther F low AC, 20-44 y")))
  
  #rp.2.2.2.1
  rp.2.2.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.2.2.1"])))
  rp.2.2.2.1.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.2.2.2.1"],prior.param2["rp.2.2.2.1"])))
  plot.rp.2.2.2.1 <- ggplot() +
    geom_histogram(data=rp.2.2.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.2.2.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nOther F high AC, 20-44 y")))
  
  #rp.2.1.1.2
  rp.2.1.1.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.1.1.2"])))
  rp.2.1.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.2.1.1.2"],prior.param2["rp.2.1.1.2"])))
  plot.rp.2.1.1.2 <- ggplot() +
    geom_histogram(data=rp.2.1.1.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.1.1.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR partner change:\nOther M low AC, 45-64 y")))
  
  #rp.2.1.2.2
  rp.2.1.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.1.1.2.2"])))
  rp.2.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.2.1.2.2"],prior.param2["rp.2.1.2.2"])))
  plot.rp.2.1.2.2 <- ggplot() +
    geom_histogram(data=rp.2.1.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.2.1.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nOther M high AC, 45-64 y")))
  
  #rp.2.2.1.2
  rp.2.2.1.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.2.1.2"])))
  rp.2.2.1.2.prior <- as.data.frame(cbind(x=seq(0,5,0.1),y=dgamma(seq(0,5,0.1),prior.param1["rp.2.2.1.2"],prior.param2["rp.2.2.1.2"])))
  plot.rp.2.2.1.2 <- ggplot() +
    geom_histogram(data=rp.2.2.1.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.2.2.1.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,5))+
    labs(x=expression(paste("RR partner change:\nOther F low AC, 45-64 y")))
  
  #rp.2.2.2.2
  rp.2.2.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.2.2.2.2"])))
  rp.2.2.2.2.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.2.2.2.2"],prior.param2["rp.2.2.2.2"])))
  plot.rp.2.2.2.2 <- ggplot() +
    geom_histogram(data=rp.2.2.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.2.2.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nOther F high AC, 45-64 y")))
  
  #rp.3.1.2.1
  rp.3.1.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.3.1.2.1"])))
  rp.3.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.3.1.2.1"],prior.param2["rp.3.1.2.1"])))
  plot.rp.3.1.2.1 <- ggplot() +
    geom_histogram(data=rp.3.1.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.3.1.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nHispanic M high AC, 20-44 y")))
  
  #rp.3.2.2.1
  rp.3.2.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.3.2.2.1"])))
  rp.3.2.2.1.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.3.2.2.1"],prior.param2["rp.3.2.2.1"])))
  plot.rp.3.2.2.1 <- ggplot() +
    geom_histogram(data=rp.3.2.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.3.2.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nHispanic F high AC, 20-44 y")))
  
  #rp.3.1.2.2
  rp.3.1.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.3.1.2.2"])))
  rp.3.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.3.1.2.2"],prior.param2["rp.3.1.2.2"])))
  plot.rp.3.1.2.2 <- ggplot() +
    geom_histogram(data=rp.3.1.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=2)+
    geom_area(data=rp.3.1.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nHispanic M high AC, 45-64 y")))
  
  #rp.3.2.2.2
  rp.3.2.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.3.2.2.2"])))
  rp.3.2.2.2.prior <- as.data.frame(cbind(x=seq(0,50,0.1),y=dgamma(seq(0,50,0.1),prior.param1["rp.3.2.2.2"],prior.param2["rp.3.2.2.2"])))
  plot.rp.3.2.2.2 <- ggplot() +
    geom_histogram(data=rp.3.2.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=1)+
    geom_area(data=rp.3.2.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,50))+
    labs(x=expression(paste("RR partner change:\nHispanic F high AC, 45-64 y")))
  
  #rp.4.1.2.1
  rp.4.1.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.4.1.2.1"])))
  rp.4.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.4.1.2.1"],prior.param2["rp.4.1.2.1"])))
  plot.rp.4.1.2.1 <- ggplot() +
    geom_histogram(data=rp.4.1.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.4.1.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nMSM high AC, 20-44 y")))
  
  #rp.4.1.2.2
  rp.4.1.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.4.1.2.2"])))
  rp.4.1.2.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.4.1.2.2"],prior.param2["rp.4.1.2.2"])))
  plot.rp.4.1.2.2 <- ggplot() +
    geom_histogram(data=rp.4.1.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.4.1.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nMSM high AC, 45-64 y")))
  
  #rp.5.1.1.1
  rp.5.1.1.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.5.1.1.1"])))
  rp.5.1.1.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.5.1.1.1"],prior.param2["rp.5.1.1.1"])))
  plot.rp.5.1.1.1 <- ggplot() +
    geom_histogram(data=rp.5.1.1.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.5.1.1.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nHIV+ MSM low AC, 20-44 y")))
  
  #rp.5.1.2.1
  rp.5.1.2.1.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.5.1.2.1"])))
  rp.5.1.2.1.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dnorm(seq(0,100,0.1),prior.param1["rp.5.1.2.1"],prior.param2["rp.5.1.2.1"])))
  plot.rp.5.1.2.1 <- ggplot() +
    geom_histogram(data=rp.5.1.2.1.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.5.1.2.1.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3) +
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,100))+
    labs(x=expression(paste("RR partner change:\nHIV+ MSM high AC, 20-44 y")))
  
  #rp.5.1.1.2
  rp.5.1.1.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.5.1.1.2"])))
  rp.5.1.1.2.prior <- as.data.frame(cbind(x=seq(0,100,0.1),y=dgamma(seq(0,100,0.1),prior.param1["rp.5.1.1.2"],prior.param2["rp.5.1.1.2"])))
  plot.rp.5.1.1.2 <- ggplot() +
    geom_histogram(data=rp.5.1.1.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=0.2)+
    geom_area(data=rp.5.1.1.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,10))+
    labs(x=expression(paste("RR partner change:\nHIV+ MSM low AC, 45-64 y")))
  
  #rp.5.1.2.2
  rp.5.1.2.2.post<-data.frame(X1=as.numeric(exp(trace.burn.thin[,"log.rp.5.1.2.2"])))
  rp.5.1.2.2.prior <- as.data.frame(cbind(x=seq(0,150,0.1),y=dnorm(seq(0,150,0.1),prior.param1["rp.5.1.2.2"],prior.param2["rp.5.1.2.2"])))
  plot.rp.5.1.2.2 <- ggplot() +
    geom_histogram(data=rp.5.1.2.2.post, aes(x=X1, y=..density..), fill="mediumpurple1",size=0.5, colour="mediumpurple1", alpha=0.6, binwidth=5)+
    geom_area(data=rp.5.1.2.2.prior,aes(x=x, y=y), fill="dimgrey", alpha=0.3)+
    theme_classic() +
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8), axis.title.x=element_text(margin=margin(t=10))) + 
    coord_cartesian(xlim=c(0,150))+
    labs(x=expression(paste("RR partner change:\nHIV+ MSM high AC, 45-64 y")))
  
  ### prepare time-varying parameters ###
  x<-apply(post.sample$theta,1,update.ctrl)
  post.screen.bc<-sapply(x,`[`,2)
  post.rep.bc<-sapply(x,`[`,1)
  post.rep.b<-unname(sapply(post.rep.bc,`[`,1))
  post.rep.c<-unname(sapply(post.rep.bc,`[`,2))
  post.scr.m1.b<-sapply(post.screen.bc,`[`,1)
  post.scr.msm1.b<-sapply(post.screen.bc,`[`,2)
  post.scr.f1.b<-sapply(post.screen.bc,`[`,3)
  post.scr.m1.c<-sapply(post.screen.bc,`[`,4)
  post.scr.msm1.c<-sapply(post.screen.bc,`[`,5)
  post.scr.f1.c<-sapply(post.screen.bc,`[`,6)
  
  #rep.symp - reporting probability
  bez.rep <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.rep.a"]), bezD=ilogit(post.sample$theta[,"logit.rep.d"]),bezB= post.rep.b, bezC=post.rep.c))
  rep.post <-reshape::melt(apply(bez.rep, 1, bezier.fun, length.out = (cal.period + 10)))
  #browser()
  rep.post$X1 <- rep.post$X1+start.year-11
  #rep.post$X1 <- rep.post$X1+start.year-11
  rep.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["rep.a"],prior.param2["rep.a"]),d=rbeta(1000,prior.param1["rep.d"],prior.param2["rep.d"]),b=runif(1000,prior.param1["rand.rep.b"],prior.param2["rand.rep.b"]),c=runif(1000,prior.param1["rand.rep.c"],prior.param2["rand.rep.c"])))
  rep.prior.bez<-matrix(apply(rep.prior.bez,1, prior.ctrl),ncol=4, byrow=TRUE)
  rep.prior <-reshape::melt(apply(rep.prior.bez, 1, bezier.fun, length.out = (cal.period + 10)))
  x<-data.frame(c(aggregate(value~X1, rep.prior, mean),aggregate(value~X1,rep.prior, min),  aggregate(value~X1, rep.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")
  
  plot.rep <- ggplot(data=rep.post) +
    geom_line(aes(x=X1, y=value, group=X2), color="thistle4")+
    geom_ribbon(data=x, aes(x=seq(start.year-10, end.year, 1),ymin=min, ymax=max), alpha=0.2)+
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) + 
    labs(x="Year", y="Reporting probability\n(identified by screening)") 
  
  #transmission RR in MSM - change in condom use/behaviour in MSM leading to changes in transmission over time
  bez.behav <- data.frame(X1=as.numeric(ilogit(post.sample$theta[,"logit.behav.lin"])))
  behav.post <-reshape::melt(apply(bez.behav, 1, behav.fun))
  behav.post$X1 <- behav.post$X1+start.year-11
  behav.prior.bez<- as.data.frame(rbeta(1000,prior.param1["behav.lin"],prior.param2["behav.lin"]))
  behav.prior <-reshape::melt(apply(behav.prior.bez, 1, behav.fun))
  x<-data.frame(c(aggregate(value~X1, behav.prior, mean),aggregate(value~X1,behav.prior, min),  aggregate(value~X1, behav.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")
  
  plot.behav <- ggplot(data=behav.post) +
    geom_line(aes(x=X1, y=value, group=X2, color=X2))+
    geom_ribbon(data=x, aes(x=(start.year-10):(end.year-1),ymin=min, ymax=max), alpha=0.2)+
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) + 
    labs(x="Year", y="Transmission\nRR in MSM")
  
  #screen.m1 - screening rate in youngest males (other)
  bez.screen.m1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.m1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.m1.d"]),bezB= post.scr.m1.b, bezC=post.scr.m1.c))
  screen.m1.post <-reshape::melt(apply(bez.screen.m1, 1, bezier.fun, length.out = (cal.period+10)))
  screen.m1.post$X1 <- screen.m1.post$X1+start.year-11
  screen.m1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.m1.a"],prior.param2["screen.m1.a"]),d=rbeta(1000,prior.param1["screen.m1.d"],prior.param2["screen.m1.d"]),b=runif(1000,prior.param1["rand.screen.m1.b"],prior.param2["rand.screen.m1.b"]),c=runif(1000,prior.param1["rand.screen.m1.c"],prior.param2["rand.screen.m1.c"])))
  screen.m1.prior.bez <- matrix(apply(screen.m1.prior.bez,1, prior.ctrl),ncol=4, byrow=TRUE)
  screen.m1.prior <-reshape::melt(apply(screen.m1.prior.bez, 1, bezier.fun, length.out = (cal.period+10)))
  x<-data.frame(c(aggregate(value~X1, screen.m1.prior, mean),aggregate(value~X1,screen.m1.prior, min),  aggregate(value~X1, screen.m1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")
  
  plot.screen.m1 <- ggplot(data=screen.m1.post) +
    geom_line(aes(x=X1, y=value, group=X2), color="steelblue1")+
    geom_ribbon(data=x, aes(x=seq(start.year-10, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) + 
    labs(title="M, other: 20-44y",x="Year", y="Screening rate") +
    coord_cartesian(ylim=c(0,1))
  
  #screen.msm1 - screening rate in youngest MSM
  bez.screen.msm1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.msm1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.msm1.d"]),bezB= post.scr.msm1.b, bezC=post.scr.msm1.c))
  screen.msm1.post <-reshape::melt(apply(bez.screen.msm1, 1, bezier.fun, length.out = (cal.period+10)))
  screen.msm1.post$X1 <- screen.msm1.post$X1+start.year-11
  screen.msm1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.msm1.a"],prior.param2["screen.msm1.a"]),d=rbeta(1000,prior.param1["screen.msm1.d"],prior.param2["screen.msm1.d"]),b=runif(1000,prior.param1["rand.screen.msm1.b"],prior.param2["rand.screen.msm1.b"]),c=runif(1000,prior.param1["rand.screen.msm1.c"],prior.param2["rand.screen.msm1.c"])))
  screen.msm1.prior.bez <- matrix(apply(screen.msm1.prior.bez,1, prior.ctrl),ncol=4, byrow=TRUE)
  screen.msm1.prior <-reshape::melt(apply(screen.msm1.prior.bez, 1, bezier.fun, length.out = (cal.period+10)))
  x<-data.frame(c(aggregate(value~X1, screen.msm1.prior, mean),aggregate(value~X1,screen.msm1.prior, min),  aggregate(value~X1, screen.msm1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")
  
  plot.screen.msm1 <- ggplot(data=screen.msm1.post) +
    geom_line(aes(x=X1, y=value, group=X2), color="steelblue1")+
    geom_ribbon(data=x, aes(x=seq(start.year-10, end.year,1),ymin=min, ymax=max), alpha=0.2)+
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) + 
    labs(title="MSM: 20-44y", x="Year", y="Screening rate") +
    coord_cartesian(ylim=c(0,1))
  
  #screen.f1 - screening rate in youngest females (other)
  bez.screen.f1 <- as.data.frame(cbind(bezA=ilogit(post.sample$theta[,"logit.screen.f1.a"]), bezD=ilogit(post.sample$theta[,"logit.screen.f1.d"]),bezB= post.scr.f1.b, bezC=post.scr.f1.c))
  screen.f1.post <-reshape::melt(apply(bez.screen.f1, 1, bezier.fun, length.out = (cal.period+10)))
  screen.f1.post$X1 <- screen.f1.post$X1+start.year-11
  screen.f1.prior.bez<- as.data.frame(cbind(a=rbeta(1000,prior.param1["screen.f1.a"],prior.param2["screen.f1.a"]),d=rbeta(1000,prior.param1["screen.f1.d"],prior.param2["screen.f1.d"]),b=runif(1000,prior.param1["rand.screen.f1.b"],prior.param2["rand.screen.f1.b"]),c=runif(1000,prior.param1["rand.screen.f1.c"],prior.param2["rand.screen.f1.c"])))
  screen.f1.prior.bez <- matrix(apply(screen.f1.prior.bez,1, prior.ctrl),ncol=4, byrow=TRUE)
  screen.f1.prior <-reshape::melt(apply(screen.f1.prior.bez, 1, bezier.fun, length.out = (cal.period+10)))
  x<-data.frame(c(aggregate(value~X1, screen.f1.prior, mean),aggregate(value~X1,screen.f1.prior, min),  aggregate(value~X1, screen.f1.prior, max)))
  x<-x[,c(1,2,4,6)]
  names(x)<-c("time", "mean", "min", "max")
  
  plot.screen.f1 <- ggplot(data=screen.f1.post) +
    geom_line(aes(x=X1, y=value, group=X2), color="steelblue1") +
    geom_ribbon(data=x, aes(x=seq(start.year-10, end.year, 1),ymin=min, ymax=max), alpha=0.2)+
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(size=8),axis.text.y=element_text(size=8), title=element_text(size=8)) + 
    labs(title="F, other: 20-44y",x="Year", y="Screening rate") +
    coord_cartesian(ylim=c(0,1))
  ### save prior/posterior plots to pdf ###
  pdf(file=file.path(output_dir, paste("syph_plots_",state, ".pdf", sep="")), width=10, height=10, paper="USr")
  grid.arrange(plot.diag.y.m, plot.diag.o.m, 
               plot.diag.y.f, plot.diag.o.f, 
               bottom=textGrob("Calibration targets", x=0.01, y=1, just="left"),
               ncol=2 )
  grid.arrange(plot.p.sec.y.m, plot.p.sec.o.m, plot.p.sec.y.f, 
               plot.p.el.y.m,  plot.p.el.o.m,  plot.p.el.y.f, 
               bottom=textGrob("Calibration targets", x=0.01, y=1, just="left"),
               ncol=3)
  grid.arrange(plot.p.msm.y, plot.p.msm.o, 
               plot.p.hiv.y, plot.p.hiv.o, 
               bottom=textGrob("Calibration targets", x=0.01, y=1, just="left"),
               ncol=2)
  grid.arrange(plot.subpop, plot.age, plot.rr.diag, plot.p.early, ncol=2)
  
  grid.arrange(plot.inc.y.m.1, plot.inc.y.m.2, plot.inc.y.m.3, plot.inc.y.msm, plot.inc.y.msmhiv, 
               plot.inc.o.m.1, plot.inc.o.m.2, plot.inc.o.m.3, plot.inc.o.msm, plot.inc.o.msmhiv, 
               plot.inc.y.f.1, plot.inc.y.f.2, plot.inc.y.f.3, blankPlot,      blankPlot,
               plot.inc.o.f.1, plot.inc.o.f.2, plot.inc.o.f.3, blankPlot,      blankPlot,   ncol=5)
  
  grid.arrange(plot.b.m,         plot.b.f,           plot.b.msm,      plot.behav,
               plot.dur.incub,   plot.dur.prim,      plot.dur.sec,    blankPlot,
               plot.dur.imm.inf, plot.dur.imm.early, plot.dur.immune, blankPlot, 
               plot.theta.1,     plot.theta.2,       plot.theta.3,    plot.theta.4, 
               plot.theta.5,     plot.theta.6,       plot.theta.7,    plot.theta.8,
               plot.epsilon.1,   plot.epsilon.2,     plot.epsilon.3,  plot.epsilon.4,
               plot.epsilon.5,   plot.theta.hiv,
               plot.pi.m,        plot.pi.f,          plot.pi.msm,     blankPlot,
               bottom=textGrob("Priors/posteriors (1/4)", x=0.01, y=1, just="left"),
               ncol=4)
  
  grid.arrange(plot.trt.prim.m,    plot.trt.prim.f,    plot.trt.prim.msm,    plot.trt.prim.msmhiv,
               plot.trt.sec.m,     plot.trt.sec.f,     plot.trt.sec.msm,     plot.trt.sec.msmhiv,
               plot.trt.lat.m,     plot.trt.lat.f,     plot.trt.lat.msm,     plot.trt.lat.msmhiv,
               plot.trt.late.m,    plot.trt.late.f,    plot.trt.late.msm,    plot.trt.late.msmhiv,
               bottom=textGrob("Priors/posteriors (2/4)", x=0.01, y=1, just="left"),
               ncol=4)
  #browser()
  grid.arrange(plot.screen.m1,     plot.screen.f1,     plot.screen.msm1,   
               plot.rr.screen.m1,  plot.rr.screen.f1,  plot.rr.screen.o.msm, 
               plot.rr.screen.m3,  plot.rr.screen.f3,  plot.rr.screen.msmhiv,
               plot.rr.screen.o.m, plot.rr.screen.o.f, plot.rr.screen.ac,            
               plot.rr.rep.symp.m, plot.rr.rep.symp.f, plot.rep,            
               bottom=textGrob("Priors/posteriors (3/4)", x=0.01, y=1, just="left"),
               ncol=3)
  
  grid.arrange(plot.rp.1.1.1.1, plot.rp.2.1.1.1, plot.c.min.m1,   plot.c.min.msm1, plot.rp.5.1.1.1, 
               plot.rp.1.1.1.2, plot.rp.2.1.1.2, plot.c.min.m2,   plot.c.min.msm2, plot.rp.5.1.1.2,
               plot.rp.1.2.1.1, plot.rp.2.2.1.1, plot.c.min.f1,   blankPlot,       blankPlot,
               plot.rp.1.2.1.2, plot.rp.2.2.1.2, plot.c.min.f2,   blankPlot,       blankPlot,
               plot.rp.1.1.2.1, plot.rp.2.1.2.1, plot.rp.3.1.2.1, plot.rp.4.1.2.1, plot.rp.5.1.2.1, 
               plot.rp.1.1.2.2, plot.rp.2.1.2.2, plot.rp.3.1.2.2, plot.rp.4.1.2.2, plot.rp.5.1.2.2, 
               plot.rp.1.2.2.1, plot.rp.2.2.2.1, plot.rp.3.2.2.1, blankPlot,       blankPlot,
               plot.rp.1.2.2.2, plot.rp.2.2.2.2, plot.rp.3.2.2.2, blankPlot,       blankPlot,
               bottom=textGrob("Priors/posteriors (4/4)", x=0.01, y=1, just="left"),
               ncol=5)
  
  if(exists("showCounterfactual") && showCounterfactual == TRUE) {
  grid.arrange(plot.diag.y.m.cf, plot.diag.o.m.cf, 
               plot.diag.y.f.cf, plot.diag.o.f.cf, 
               bottom=textGrob("Calibration targets", x=0.01, y=1, just="left"),
               ncol=2 )
  grid.arrange(plot.inc.y.m.1.cf, plot.inc.y.m.2.cf, plot.inc.y.m.3.cf, plot.inc.y.msm.cf, plot.inc.y.msmhiv.cf, 
               plot.inc.o.m.1.cf, plot.inc.o.m.2.cf, plot.inc.o.m.3.cf, plot.inc.o.msm.cf, plot.inc.o.msmhiv.cf, 
               plot.inc.y.f.1.cf, plot.inc.y.f.2.cf, plot.inc.y.f.3.cf, blankPlot,      blankPlot,
               plot.inc.o.f.1.cf, plot.inc.o.f.2.cf, plot.inc.o.f.3.cf, blankPlot,      blankPlot,   ncol=5)
  }
  dev.off()
  
}
