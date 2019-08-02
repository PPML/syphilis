###################################################
###   function to calculate prior likelihoods   ###
###################################################
#reads in current value for each parameter (contained in theta), and calculates likelihood
#asssumed distribution type for each parameter, and associated parameters describing each distribution contained in input file 'priors
#prior.param1 = 1st parameter for distribution; prior.param2=2nd param for distribution 
#if add or remove parameters from calibration, need to modify this code

old_prior <- function(theta) {
  dbeta(ilogit(theta["logit.epsilon.1"]),prior.param1["epsilon.1"],prior.param2["epsilon.1"])* 
  dbeta(ilogit(theta["logit.epsilon.2"]),prior.param1["epsilon.2"],prior.param2["epsilon.2"])*
  dbeta(ilogit(theta["logit.epsilon.3"]),prior.param1["epsilon.3"],prior.param2["epsilon.3"])*
  dbeta(ilogit(theta["logit.epsilon.4"]),prior.param1["epsilon.4"],prior.param2["epsilon.4"])*
  dbeta(ilogit(theta["logit.epsilon.5"]),prior.param1["epsilon.5"],prior.param2["epsilon.5"])*
  dgamma(exp(theta["log.rr.screen.m1"]),prior.param1["rr.screen.m1"],prior.param2["rr.screen.m1"]) *
  dgamma(exp(theta["log.rr.screen.m3"]),prior.param1["rr.screen.m3"],prior.param2["rr.screen.m3"]) *
  dgamma(exp(theta["log.rr.screen.o.m"]),prior.param1["rr.screen.o.m"],prior.param2["rr.screen.o.m"]) *
  dgamma(exp(theta["log.rr.screen.msmhiv"]),prior.param1["rr.screen.msmhiv"],prior.param2["rr.screen.msmhiv"]) *
  dgamma(exp(theta["log.rr.screen.o.msm"]),prior.param1["rr.screen.o.msm"],prior.param2["rr.screen.o.msm"]) *
  dgamma(exp(theta["log.rr.screen.f1"]),prior.param1["rr.screen.f1"],prior.param2["rr.screen.f1"]) *
  dgamma(exp(theta["log.rr.screen.f3"]),prior.param1["rr.screen.f3"],prior.param2["rr.screen.f3"]) *
  dgamma(exp(theta["log.rr.screen.o.f"]),prior.param1["rr.screen.o.f"],prior.param2["rr.screen.o.f"]) *
  dgamma(exp(theta["log.rr.screen.ac"]),prior.param1["rr.screen.ac"],prior.param2["rr.screen.ac"]) *
  dgamma(exp(theta["log.c.min.m.1"]),prior.param1["c.min.m1"],prior.param2["c.min.m1"]) * 
  dgamma(exp(theta["log.c.min.m.2"]),prior.param1["c.min.m2"],prior.param2["c.min.m2"]) *
  dgamma(exp(theta["log.c.min.f.1"]),prior.param1["c.min.f1"],prior.param2["c.min.f1"]) * 
  dgamma(exp(theta["log.c.min.f.2"]),prior.param1["c.min.f2"],prior.param2["c.min.f2"])*
  dgamma(exp(theta["log.c.min.msm.1"]),prior.param1["c.min.msm1"],prior.param2["c.min.msm1"]) * 
  dgamma(exp(theta["log.c.min.msm.2"]),prior.param1["c.min.msm2"],prior.param2["c.min.msm2"]) *
  dbeta(ilogit(theta["logit.b.m"]), prior.param1["b.m"],prior.param2["b.m"]) *
  dbeta(ilogit(theta["logit.b.f"]), prior.param1["b.f"],prior.param2["b.f"]) *
  dbeta(ilogit(theta["logit.b.msm"]), prior.param1["b.msm"],prior.param2["b.msm"])*
  dnorm(exp(theta["log.dur.incub"]),prior.param1["dur.incub"],prior.param2["dur.incub"]) * 
  dnorm(exp(theta["log.dur.prim"]),prior.param1["dur.prim"],prior.param2["dur.prim"]) * 
  dnorm(exp(theta["log.dur.sec"]),prior.param1["dur.sec"],prior.param2["dur.sec"]) *
  dnorm(exp(theta["log.dur.imm.inf"]),prior.param1["dur.imm.inf"],prior.param2["dur.imm.inf"]) * 
  dunif(exp(theta["log.dur.imm.early"]),prior.param1["dur.imm.early"],prior.param2["dur.imm.early"]) * 
  dnorm(exp(theta["log.dur.immune"]),prior.param1["dur.immune"],prior.param2["dur.immune"]) *
  dbeta(ilogit(theta["logit.p.trt.prim.m"]),prior.param1["p.trt.prim.m"],prior.param2["p.trt.prim.m"]) * 
  dbeta(ilogit(theta["logit.p.trt.prim.f"]),prior.param1["p.trt.prim.f"],prior.param2["p.trt.prim.f"]) * 
  dbeta(ilogit(theta["logit.p.trt.prim.msm"]),prior.param1["p.trt.prim.msm"],prior.param2["p.trt.prim.msm"]) *
  dbeta(ilogit(theta["logit.p.trt.prim.msmhiv"]),prior.param1["p.trt.prim.msmhiv"],prior.param2["p.trt.prim.msmhiv"]) *
  dbeta(ilogit(theta["logit.p.trt.sec.m"]),prior.param1["p.trt.sec.m"],prior.param2["p.trt.sec.m"]) * 
  dbeta(ilogit(theta["logit.p.trt.sec.f"]),prior.param1["p.trt.sec.f"],prior.param2["p.trt.sec.f"]) * 
  dbeta(ilogit(theta["logit.p.trt.sec.msm"]),prior.param1["p.trt.sec.msm"],prior.param2["p.trt.sec.msm"]) *
  dbeta(ilogit(theta["logit.p.trt.sec.msmhiv"]),prior.param1["p.trt.sec.msmhiv"],prior.param2["p.trt.sec.msmhiv"]) *
  dbeta(ilogit(theta["logit.p.trt.lat.m"]),prior.param1["p.trt.lat.m"],prior.param2["p.trt.lat.m"]) * 
  dbeta(ilogit(theta["logit.p.trt.lat.f"]),prior.param1["p.trt.lat.f"],prior.param2["p.trt.lat.f"]) * 
  dbeta(ilogit(theta["logit.p.trt.lat.msm"]),prior.param1["p.trt.lat.msm"],prior.param2["p.trt.lat.msm"]) *
  dbeta(ilogit(theta["logit.p.trt.lat.msmhiv"]),prior.param1["p.trt.lat.msmhiv"],prior.param2["p.trt.lat.msmhiv"]) *
  dbeta(ilogit(theta["logit.p.trt.late.m"]),prior.param1["p.trt.late.m"],prior.param2["p.trt.late.m"]) * 
  dbeta(ilogit(theta["logit.p.trt.late.f"]),prior.param1["p.trt.late.f"],prior.param2["p.trt.late.f"]) * 
  dbeta(ilogit(theta["logit.p.trt.late.msm"]),prior.param1["p.trt.late.msm"],prior.param2["p.trt.late.msm"]) *
  dbeta(ilogit(theta["logit.p.trt.late.msmhiv"]),prior.param1["p.trt.late.msmhiv"],prior.param2["p.trt.late.msmhiv"]) *
  dbeta(ilogit(theta["logit.theta.1"]),prior.param1["theta.1"],prior.param2["theta.1"]) * 
  dbeta(ilogit(theta["logit.theta.2"]),prior.param1["theta.2"],prior.param2["theta.2"]) * 
  dbeta(ilogit(theta["logit.theta.3"]),prior.param1["theta.3"],prior.param2["theta.3"]) * 
  dbeta(ilogit(theta["logit.theta.4"]),prior.param1["theta.4"],prior.param2["theta.4"])*
  dbeta(ilogit(theta["logit.theta.5"]),prior.param1["theta.5"],prior.param2["theta.5"]) * 
  dbeta(ilogit(theta["logit.theta.6"]),prior.param1["theta.6"],prior.param2["theta.6"]) * 
  dbeta(ilogit(theta["logit.theta.7"]),prior.param1["theta.7"],prior.param2["theta.7"]) * 
  dbeta(ilogit(theta["logit.theta.8"]),prior.param1["theta.8"],prior.param2["theta.8"]) *
  dbeta(ilogit(theta["logit.theta.hiv"]),prior.param1["theta.hiv"],prior.param2["theta.hiv"]) *
  dbeta(ilogit(theta["logit.pi.m"]),prior.param1["pi.m"],prior.param2["pi.m"]) * 
  dbeta(ilogit(theta["logit.pi.f"]),prior.param1["pi.f"],prior.param2["pi.f"]) * 
  dbeta(ilogit(theta["logit.pi.msm"]),prior.param1["pi.msm"],prior.param2["pi.msm"]) *
  dgamma(exp(theta["log.rp.1.1.1.1"]),prior.param1["rp.1.1.1.1"], prior.param2["rp.1.1.1.1"])*
  dnorm(exp(theta["log.rp.1.1.2.1"]),prior.param1["rp.1.1.2.1"], prior.param2["rp.1.1.2.1"])*
  dgamma(exp(theta["log.rp.1.2.1.1"]),prior.param1["rp.1.2.1.1"], prior.param2["rp.1.2.1.1"])*
  dgamma(exp(theta["log.rp.1.2.2.1"]),prior.param1["rp.1.2.2.1"], prior.param2["rp.1.2.2.1"])*
  dgamma(exp(theta["log.rp.1.1.1.2"]),prior.param1["rp.1.1.1.2"], prior.param2["rp.1.1.1.2"])*
  dnorm(exp(theta["log.rp.1.1.2.2"]),prior.param1["rp.1.1.2.2"], prior.param2["rp.1.1.2.2"])*
  dgamma(exp(theta["log.rp.1.2.1.2"]),prior.param1["rp.1.2.1.2"], prior.param2["rp.1.2.1.2"])* 
  dgamma(exp(theta["log.rp.1.2.2.2"]),prior.param1["rp.1.2.2.2"], prior.param2["rp.1.2.2.2"])*
  dgamma(exp(theta["log.rp.2.1.1.1"]),prior.param1["rp.2.1.1.1"], prior.param2["rp.2.1.1.1"])*
  dnorm(exp(theta["log.rp.2.1.2.1"]),prior.param1["rp.2.1.2.1"], prior.param2["rp.2.1.2.1"])*
  dgamma(exp(theta["log.rp.2.2.1.1"]),prior.param1["rp.2.2.1.1"], prior.param2["rp.2.2.1.1"])*
  dgamma(exp(theta["log.rp.2.2.2.1"]),prior.param1["rp.2.2.2.1"], prior.param2["rp.2.2.2.1"])*
  dgamma(exp(theta["log.rp.2.1.1.2"]),prior.param1["rp.2.1.1.2"], prior.param2["rp.2.1.1.2"])* 
  dnorm(exp(theta["log.rp.2.1.2.2"]),prior.param1["rp.2.1.2.2"], prior.param2["rp.2.1.2.2"])*
  dgamma(exp(theta["log.rp.2.2.1.2"]),prior.param1["rp.2.2.1.2"], prior.param2["rp.2.2.1.2"])* 
  dgamma(exp(theta["log.rp.2.2.2.2"]),prior.param1["rp.2.2.2.2"], prior.param2["rp.2.2.2.2"])*
  dgamma(exp(theta["log.rp.3.1.2.1"]),prior.param1["rp.3.1.2.1"], prior.param2["rp.3.1.2.1"])*
  dgamma(exp(theta["log.rp.3.2.2.1"]),prior.param1["rp.3.2.2.1"], prior.param2["rp.3.2.2.1"])*
  dgamma(exp(theta["log.rp.3.1.2.2"]),prior.param1["rp.3.1.2.2"], prior.param2["rp.3.1.2.2"])* 
  dgamma(exp(theta["log.rp.3.2.2.2"]),prior.param1["rp.3.2.2.2"], prior.param2["rp.3.2.2.2"])*
  dnorm(exp(theta["log.rp.4.1.2.1"]),prior.param1["rp.4.1.2.1"], prior.param2["rp.4.1.2.1"])*
  dnorm(exp(theta["log.rp.4.1.2.2"]),prior.param1["rp.4.1.2.2"], prior.param2["rp.4.1.2.2"])*
  dgamma(exp(theta["log.rp.5.1.1.1"]),prior.param1["rp.5.1.1.1"], prior.param2["rp.5.1.1.1"])*
  dnorm(exp(theta["log.rp.5.1.2.1"]),prior.param1["rp.5.1.2.1"], prior.param2["rp.5.1.2.1"])*
  dgamma(exp(theta["log.rp.5.1.1.2"]),prior.param1["rp.5.1.1.2"], prior.param2["rp.5.1.1.2"])*
  dnorm(exp(theta["log.rp.5.1.2.2"]),prior.param1["rp.5.1.2.2"], prior.param2["rp.5.1.2.2"])*
  dbeta(ilogit(theta["logit.rep.a"]),prior.param1["rep.a"],prior.param2["rep.a"]) * 
  dbeta(ilogit(theta["logit.rep.d"]),prior.param1["rep.d"],prior.param2["rep.d"])*
  dunif(ilogit(theta["logit.rand.rep.b"]),prior.param1["rand.rep.b"],prior.param2["rand.rep.b"]) * 
  dunif(ilogit(theta["logit.rand.rep.c"]),prior.param1["rand.rep.c"],prior.param2["rand.rep.c"])*
  dbeta(ilogit(theta["logit.rr.rep.symp.m"]),prior.param1["rr.rep.symp.m"],prior.param2["rr.rep.symp.m"])*
  dbeta(ilogit(theta["logit.rr.rep.symp.f"]),prior.param1["rr.rep.symp.f"],prior.param2["rr.rep.symp.f"])*
  dunif(ilogit(theta["logit.behav.lin"]),prior.param1["behav.lin"],prior.param2["behav.lin"])*
  dbeta(ilogit(theta["logit.screen.m1.a"]),prior.param1["screen.m1.a"],prior.param2["screen.m1.a"]) * 
  dbeta(ilogit(theta["logit.screen.m1.d"]),prior.param1["screen.m1.d"],prior.param2["screen.m1.d"])*
  dunif(ilogit(theta["logit.rand.screen.m1.b"]),prior.param1["rand.screen.m1.b"],prior.param2["rand.screen.m1.b"]) * 
  dunif(ilogit(theta["logit.rand.screen.m1.c"]),prior.param1["rand.screen.m1.c"],prior.param2["rand.screen.m1.c"])*
  dbeta(ilogit(theta["logit.screen.msm1.a"]),prior.param1["screen.msm1.a"],prior.param2["screen.msm1.a"]) * 
  dbeta(ilogit(theta["logit.screen.msm1.d"]),prior.param1["screen.msm1.d"],prior.param2["screen.msm1.d"])*
  dunif(ilogit(theta["logit.rand.screen.msm1.b"]),prior.param1["rand.screen.msm1.b"],prior.param2["rand.screen.msm1.b"]) * 
  dunif(ilogit(theta["logit.rand.screen.msm1.c"]),prior.param1["rand.screen.msm1.c"],prior.param2["rand.screen.msm1.c"])*
  dbeta(ilogit(theta["logit.screen.f1.a"]),prior.param1["screen.f1.a"],prior.param2["screen.f1.a"]) * 
  dbeta(ilogit(theta["logit.screen.f1.d"]),prior.param1["screen.f1.d"],prior.param2["screen.f1.d"])*
  dunif(ilogit(theta["logit.rand.screen.f1.b"]),prior.param1["rand.screen.f1.b"],prior.param2["rand.screen.f1.b"]) * 
  dunif(ilogit(theta["logit.rand.screen.f1.c"]),prior.param1["rand.screen.f1.c"],prior.param2["rand.screen.f1.c"])
}



prior_components <- function(theta) {
return(c(
  dbeta(ilogit(theta["logit.epsilon.1"]),prior.param1["epsilon.1"],prior.param2["epsilon.1"], log=T),
  dbeta(ilogit(theta["logit.epsilon.2"]),prior.param1["epsilon.2"],prior.param2["epsilon.2"], log=T),
  dbeta(ilogit(theta["logit.epsilon.3"]),prior.param1["epsilon.3"],prior.param2["epsilon.3"], log=T),
  dbeta(ilogit(theta["logit.epsilon.4"]),prior.param1["epsilon.4"],prior.param2["epsilon.4"], log=T),
  dbeta(ilogit(theta["logit.epsilon.5"]),prior.param1["epsilon.5"],prior.param2["epsilon.5"], log=T),
  dgamma(exp(theta["log.rr.screen.m1"]),prior.param1["rr.screen.m1"],prior.param2["rr.screen.m1"], log=T),
  dgamma(exp(theta["log.rr.screen.m3"]),prior.param1["rr.screen.m3"],prior.param2["rr.screen.m3"], log=T),
  dgamma(exp(theta["log.rr.screen.o.m"]),prior.param1["rr.screen.o.m"],prior.param2["rr.screen.o.m"], log=T),
  dgamma(exp(theta["log.rr.screen.msmhiv"]),prior.param1["rr.screen.msmhiv"],prior.param2["rr.screen.msmhiv"], log=T),
  dgamma(exp(theta["log.rr.screen.o.msm"]),prior.param1["rr.screen.o.msm"],prior.param2["rr.screen.o.msm"], log=T),
  dgamma(exp(theta["log.rr.screen.f1"]),prior.param1["rr.screen.f1"],prior.param2["rr.screen.f1"], log=T),
  dgamma(exp(theta["log.rr.screen.f3"]),prior.param1["rr.screen.f3"],prior.param2["rr.screen.f3"], log=T),
  dgamma(exp(theta["log.rr.screen.o.f"]),prior.param1["rr.screen.o.f"],prior.param2["rr.screen.o.f"], log=T),
  dgamma(exp(theta["log.rr.screen.ac"]),prior.param1["rr.screen.ac"],prior.param2["rr.screen.ac"], log=T),
  dgamma(exp(theta["log.c.min.m.1"]),prior.param1["c.min.m1"],prior.param2["c.min.m1"], log=T),
  dgamma(exp(theta["log.c.min.m.2"]),prior.param1["c.min.m2"],prior.param2["c.min.m2"], log=T),
  dgamma(exp(theta["log.c.min.f.1"]),prior.param1["c.min.f1"],prior.param2["c.min.f1"], log=T),
  dgamma(exp(theta["log.c.min.f.2"]),prior.param1["c.min.f2"],prior.param2["c.min.f2"], log=T),
  dgamma(exp(theta["log.c.min.msm.1"]),prior.param1["c.min.msm1"],prior.param2["c.min.msm1"], log=T),
  dgamma(exp(theta["log.c.min.msm.2"]),prior.param1["c.min.msm2"],prior.param2["c.min.msm2"], log=T),
  dbeta(ilogit(theta["logit.b.m"]), prior.param1["b.m"],prior.param2["b.m"], log=T),
  dbeta(ilogit(theta["logit.b.f"]), prior.param1["b.f"],prior.param2["b.f"], log=T),
  dbeta(ilogit(theta["logit.b.msm"]), prior.param1["b.msm"],prior.param2["b.msm"], log=T),
  dnorm(exp(theta["log.dur.incub"]),prior.param1["dur.incub"],prior.param2["dur.incub"], log=T),
  dnorm(exp(theta["log.dur.prim"]),prior.param1["dur.prim"],prior.param2["dur.prim"], log=T),
  dnorm(exp(theta["log.dur.sec"]),prior.param1["dur.sec"],prior.param2["dur.sec"], log=T),
  dnorm(exp(theta["log.dur.imm.inf"]),prior.param1["dur.imm.inf"],prior.param2["dur.imm.inf"], log=T),
  dunif(exp(theta["log.dur.imm.early"]),prior.param1["dur.imm.early"],prior.param2["dur.imm.early"], log=T),
  dnorm(exp(theta["log.dur.immune"]),prior.param1["dur.immune"],prior.param2["dur.immune"], log=T),
  dbeta(ilogit(theta["logit.p.trt.prim.m"]),prior.param1["p.trt.prim.m"],prior.param2["p.trt.prim.m"], log=T),
  dbeta(ilogit(theta["logit.p.trt.prim.f"]),prior.param1["p.trt.prim.f"],prior.param2["p.trt.prim.f"], log=T),
  dbeta(ilogit(theta["logit.p.trt.prim.msm"]),prior.param1["p.trt.prim.msm"],prior.param2["p.trt.prim.msm"], log=T),
  dbeta(ilogit(theta["logit.p.trt.prim.msmhiv"]),prior.param1["p.trt.prim.msmhiv"],prior.param2["p.trt.prim.msmhiv"], log=T),
  dbeta(ilogit(theta["logit.p.trt.sec.m"]),prior.param1["p.trt.sec.m"],prior.param2["p.trt.sec.m"], log=T),
  dbeta(ilogit(theta["logit.p.trt.sec.f"]),prior.param1["p.trt.sec.f"],prior.param2["p.trt.sec.f"], log=T),
  dbeta(ilogit(theta["logit.p.trt.sec.msm"]),prior.param1["p.trt.sec.msm"],prior.param2["p.trt.sec.msm"], log=T),
  dbeta(ilogit(theta["logit.p.trt.sec.msmhiv"]),prior.param1["p.trt.sec.msmhiv"],prior.param2["p.trt.sec.msmhiv"], log=T),
  dbeta(ilogit(theta["logit.p.trt.lat.m"]),prior.param1["p.trt.lat.m"],prior.param2["p.trt.lat.m"], log=T),
  dbeta(ilogit(theta["logit.p.trt.lat.f"]),prior.param1["p.trt.lat.f"],prior.param2["p.trt.lat.f"], log=T),
  dbeta(ilogit(theta["logit.p.trt.lat.msm"]),prior.param1["p.trt.lat.msm"],prior.param2["p.trt.lat.msm"], log=T),
  dbeta(ilogit(theta["logit.p.trt.lat.msmhiv"]),prior.param1["p.trt.lat.msmhiv"],prior.param2["p.trt.lat.msmhiv"], log=T),
  dbeta(ilogit(theta["logit.p.trt.late.m"]),prior.param1["p.trt.late.m"],prior.param2["p.trt.late.m"], log=T),
  dbeta(ilogit(theta["logit.p.trt.late.f"]),prior.param1["p.trt.late.f"],prior.param2["p.trt.late.f"], log=T),
  dbeta(ilogit(theta["logit.p.trt.late.msm"]),prior.param1["p.trt.late.msm"],prior.param2["p.trt.late.msm"], log=T),
  dbeta(ilogit(theta["logit.p.trt.late.msmhiv"]),prior.param1["p.trt.late.msmhiv"],prior.param2["p.trt.late.msmhiv"], log=T),
  dbeta(ilogit(theta["logit.theta.1"]),prior.param1["theta.1"],prior.param2["theta.1"], log=T),
  dbeta(ilogit(theta["logit.theta.2"]),prior.param1["theta.2"],prior.param2["theta.2"], log=T),
  dbeta(ilogit(theta["logit.theta.3"]),prior.param1["theta.3"],prior.param2["theta.3"], log=T),
  dbeta(ilogit(theta["logit.theta.4"]),prior.param1["theta.4"],prior.param2["theta.4"], log=T),
  dbeta(ilogit(theta["logit.theta.5"]),prior.param1["theta.5"],prior.param2["theta.5"], log=T),
  dbeta(ilogit(theta["logit.theta.6"]),prior.param1["theta.6"],prior.param2["theta.6"], log=T),
  dbeta(ilogit(theta["logit.theta.7"]),prior.param1["theta.7"],prior.param2["theta.7"], log=T),
  dbeta(ilogit(theta["logit.theta.8"]),prior.param1["theta.8"],prior.param2["theta.8"], log=T),
  dbeta(ilogit(theta["logit.theta.hiv"]),prior.param1["theta.hiv"],prior.param2["theta.hiv"], log=T),
  dbeta(ilogit(theta["logit.pi.m"]),prior.param1["pi.m"],prior.param2["pi.m"], log=T),
  dbeta(ilogit(theta["logit.pi.f"]),prior.param1["pi.f"],prior.param2["pi.f"], log=T),
  dbeta(ilogit(theta["logit.pi.msm"]),prior.param1["pi.msm"],prior.param2["pi.msm"], log=T),
  dgamma(exp(theta["log.rp.1.1.1.1"]),prior.param1["rp.1.1.1.1"], prior.param2["rp.1.1.1.1"], log=T),
  dnorm(exp(theta["log.rp.1.1.2.1"]),prior.param1["rp.1.1.2.1"], prior.param2["rp.1.1.2.1"], log=T),
  dgamma(exp(theta["log.rp.1.2.1.1"]),prior.param1["rp.1.2.1.1"], prior.param2["rp.1.2.1.1"], log=T),
  dgamma(exp(theta["log.rp.1.2.2.1"]),prior.param1["rp.1.2.2.1"], prior.param2["rp.1.2.2.1"], log=T),
  dgamma(exp(theta["log.rp.1.1.1.2"]),prior.param1["rp.1.1.1.2"], prior.param2["rp.1.1.1.2"], log=T),
  dnorm(exp(theta["log.rp.1.1.2.2"]),prior.param1["rp.1.1.2.2"], prior.param2["rp.1.1.2.2"], log=T),
  dgamma(exp(theta["log.rp.1.2.1.2"]),prior.param1["rp.1.2.1.2"], prior.param2["rp.1.2.1.2"], log=T),
  dgamma(exp(theta["log.rp.1.2.2.2"]),prior.param1["rp.1.2.2.2"], prior.param2["rp.1.2.2.2"], log=T),
  dgamma(exp(theta["log.rp.2.1.1.1"]),prior.param1["rp.2.1.1.1"], prior.param2["rp.2.1.1.1"], log=T),
  dnorm(exp(theta["log.rp.2.1.2.1"]),prior.param1["rp.2.1.2.1"], prior.param2["rp.2.1.2.1"], log=T),
  dgamma(exp(theta["log.rp.2.2.1.1"]),prior.param1["rp.2.2.1.1"], prior.param2["rp.2.2.1.1"], log=T),
  dgamma(exp(theta["log.rp.2.2.2.1"]),prior.param1["rp.2.2.2.1"], prior.param2["rp.2.2.2.1"], log=T),
  dgamma(exp(theta["log.rp.2.1.1.2"]),prior.param1["rp.2.1.1.2"], prior.param2["rp.2.1.1.2"], log=T),
  dnorm(exp(theta["log.rp.2.1.2.2"]),prior.param1["rp.2.1.2.2"], prior.param2["rp.2.1.2.2"], log=T),
  dgamma(exp(theta["log.rp.2.2.1.2"]),prior.param1["rp.2.2.1.2"], prior.param2["rp.2.2.1.2"], log=T),
  dgamma(exp(theta["log.rp.2.2.2.2"]),prior.param1["rp.2.2.2.2"], prior.param2["rp.2.2.2.2"], log=T),
  dgamma(exp(theta["log.rp.3.1.2.1"]),prior.param1["rp.3.1.2.1"], prior.param2["rp.3.1.2.1"], log=T),
  dgamma(exp(theta["log.rp.3.2.2.1"]),prior.param1["rp.3.2.2.1"], prior.param2["rp.3.2.2.1"], log=T),
  dgamma(exp(theta["log.rp.3.1.2.2"]),prior.param1["rp.3.1.2.2"], prior.param2["rp.3.1.2.2"], log=T),
  dgamma(exp(theta["log.rp.3.2.2.2"]),prior.param1["rp.3.2.2.2"], prior.param2["rp.3.2.2.2"], log=T),
  dnorm(exp(theta["log.rp.4.1.2.1"]),prior.param1["rp.4.1.2.1"], prior.param2["rp.4.1.2.1"], log=T),
  dnorm(exp(theta["log.rp.4.1.2.2"]),prior.param1["rp.4.1.2.2"], prior.param2["rp.4.1.2.2"], log=T),
  dgamma(exp(theta["log.rp.5.1.1.1"]),prior.param1["rp.5.1.1.1"], prior.param2["rp.5.1.1.1"], log=T),
  dnorm(exp(theta["log.rp.5.1.2.1"]),prior.param1["rp.5.1.2.1"], prior.param2["rp.5.1.2.1"], log=T),
  dgamma(exp(theta["log.rp.5.1.1.2"]),prior.param1["rp.5.1.1.2"], prior.param2["rp.5.1.1.2"], log=T),
  dnorm(exp(theta["log.rp.5.1.2.2"]),prior.param1["rp.5.1.2.2"], prior.param2["rp.5.1.2.2"], log=T),
  dbeta(ilogit(theta["logit.rep.a"]),prior.param1["rep.a"],prior.param2["rep.a"], log=T),
  dbeta(ilogit(theta["logit.rep.d"]),prior.param1["rep.d"],prior.param2["rep.d"], log=T),
  dunif(ilogit(theta["logit.rand.rep.b"]),prior.param1["rand.rep.b"],prior.param2["rand.rep.b"], log=T),
  dunif(ilogit(theta["logit.rand.rep.c"]),prior.param1["rand.rep.c"],prior.param2["rand.rep.c"], log=T),
  dbeta(ilogit(theta["logit.rr.rep.symp.m"]),prior.param1["rr.rep.symp.m"],prior.param2["rr.rep.symp.m"], log=T),
  dbeta(ilogit(theta["logit.rr.rep.symp.f"]),prior.param1["rr.rep.symp.f"],prior.param2["rr.rep.symp.f"], log=T),
  dunif(ilogit(theta["logit.behav.lin"]),prior.param1["behav.lin"],prior.param2["behav.lin"], log=T),
  dbeta(ilogit(theta["logit.screen.m1.a"]),prior.param1["screen.m1.a"],prior.param2["screen.m1.a"], log=T),
  dbeta(ilogit(theta["logit.screen.m1.d"]),prior.param1["screen.m1.d"],prior.param2["screen.m1.d"], log=T),
  dunif(ilogit(theta["logit.rand.screen.m1.b"]),prior.param1["rand.screen.m1.b"],prior.param2["rand.screen.m1.b"], log=T),
  dunif(ilogit(theta["logit.rand.screen.m1.c"]),prior.param1["rand.screen.m1.c"],prior.param2["rand.screen.m1.c"], log=T),
  dbeta(ilogit(theta["logit.screen.msm1.a"]),prior.param1["screen.msm1.a"],prior.param2["screen.msm1.a"], log=T),
  dbeta(ilogit(theta["logit.screen.msm1.d"]),prior.param1["screen.msm1.d"],prior.param2["screen.msm1.d"], log=T),
  dunif(ilogit(theta["logit.rand.screen.msm1.b"]),prior.param1["rand.screen.msm1.b"],prior.param2["rand.screen.msm1.b"], log=T),
  dunif(ilogit(theta["logit.rand.screen.msm1.c"]),prior.param1["rand.screen.msm1.c"],prior.param2["rand.screen.msm1.c"], log=T),
  dbeta(ilogit(theta["logit.screen.f1.a"]),prior.param1["screen.f1.a"],prior.param2["screen.f1.a"], log=T),
  dbeta(ilogit(theta["logit.screen.f1.d"]),prior.param1["screen.f1.d"],prior.param2["screen.f1.d"], log=T),
  dunif(ilogit(theta["logit.rand.screen.f1.b"]),prior.param1["rand.screen.f1.b"],prior.param2["rand.screen.f1.b"], log=T),
  dunif(ilogit(theta["logit.rand.screen.f1.c"]),prior.param1["rand.screen.f1.c"],prior.param2["rand.screen.f1.c"], log=T)
))
}
