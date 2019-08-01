##################################################################
### function to use model outputs to calculate log likelihoods ###
##################################################################
model.epi.loglik <- function(theta) {
  pred <- prediction.epi(theta) #run the model and produce required outputs
  #calculate likelihoods
  #for beta parameter estimates, using model estimate as mean, data estimates for variance
  beta.params.age <- estBetaParamsVar(pred[["age.assort"]][c(1,3)],age.dist.dat$var)  #beta distn parameter estimates for age assortativity
  beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM
  beta.params.hiv <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.hiv"]])), p.hiv.cases.var) #beta distn parameter estimates for proportion of MSM cases who are HIV+
  beta.params.sec <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.sec"]])), p.sec.var) #beta distn parameter estimates for proportion of early cases diagnosed with secondary
  beta.params.el <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.el"]])), p.el.var) #beta distn parameter estimates for proportion of early cases diagnosed with early latent
  beta.params.early <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.early"]])), p.early.dat$var) #beta distn parameter estimates for proportion of ALL cases (early and late latent) that are early (P,S or early latent)
  
  #log likelihoods
  ll.age <- sum(dbeta(x=age.dist.dat$p.same.age,beta.params.age$alpha, beta.params.age$beta, log=TRUE))/2   #fitting to overall pop estimates only for young males and females
  ll.diag.m <- sum(dnorm(x=diag.age.sex.rate.m, mean=as.numeric(unlist(pred[["prev"]][["early.inf.rate.m"]])), sd=diag.age.sex.rate.m.sd, log=TRUE)) #log likelihood for reported male early syphilis rates
  ll.diag.f.y <- sum(dnorm(x=diag.age.sex.rate.f.y, mean=as.numeric(unlist(pred[["prev"]][["early.inf.rate.f.y"]])), sd=diag.age.sex.rate.f.y.sd, log=TRUE)) #log likelihood for reported young female early syphilis rates
  ll.diag.f.o <- sum(dnorm(x=diag.age.sex.rate.f.o, mean=as.numeric(unlist(pred[["prev"]][["early.inf.rate.f.o"]])), sd=diag.age.sex.rate.f.o.sd, log=TRUE))/4  #log likelihood for reported old female early syphilis rates
  ll.diag.subpop.m  <- sum(dnorm(x=rr.diag.subpop[1:4], mean=as.numeric(unlist(pred[["prev"]][["diag.rr"]]))[1:4],sd=rr.diag.subpop.sd[1:4], log=TRUE))/5 #log likelihood for reported case RR (black, Hispanic) males
  ll.diag.subpop.f  <- sum(dnorm(x=rr.diag.subpop[5:8], mean=as.numeric(unlist(pred[["prev"]][["diag.rr"]]))[5:8],sd=rr.diag.subpop.sd[5:8], log=TRUE))/2 #log likelihood for reported case RR (black, Hispanic) females
  ll.msm <- sum(dbeta(x=p.msm.cases,beta.params.msm$alpha, beta.params.msm$beta, log=TRUE))/2 #log likelihood for proportion of male cases in MSM
  ll.hiv <- sum(dbeta(x=p.hiv.cases,beta.params.hiv$alpha, beta.params.hiv$beta, log=TRUE))/3 #log likelihood for proportion of MSM cases with HIV coinfection
  ll.sec <- sum(dbeta(x=p.sec,beta.params.sec$alpha, beta.params.sec$beta, log=TRUE))/50 #log likelihood for proportion of early syphilis cases with secondary infection
  ll.el <- sum(dbeta(x=p.el,beta.params.el$alpha, beta.params.el$beta, log=TRUE))/50 #log likelihood for proportion of early syphilis cases with early latent infection
  ll.early <- sum(dbeta(x=p.early.dat$mean,beta.params.early$alpha, beta.params.early$beta, log=TRUE)) #log likelihood for proportion of all syphilis cases with early infection (primary, secondary or early latent)

  if (exists("overweight_early_lik")) { ll.early <- ll.early * 100 }
  if (exists("overweight_msm_lik")) { ll.msm <- ll.msm * 100 } 

  ll<-sum(ll.diag.m, ll.diag.f.y, ll.diag.f.o, ll.msm, ll.hiv, ll.el, ll.sec, ll.diag.subpop.f, ll.diag.subpop.m, ll.diag.subpop.f, ll.early)
  ll[is.na(ll)]<-(-1e20)

  
	c(unlist(pred), ll.age=ll.age, ll.diag.m=ll.diag.m, ll.diag.f.y=ll.diag.f.y, ll.diag.f.o=ll.diag.f.o,ll.rr.subpop.m= ll.diag.subpop.m,ll.diag.subpop.f=ll.diag.subpop.f ,ll.msm=ll.msm,ll.hiv=ll.hiv,ll.sec=ll.sec, ll.el=ll.el,ll.early = ll.early, ll=ll)
  
}

