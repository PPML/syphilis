devtools::load_all()
library(ggplot2)

state <- 'MA'
load.start()

# pred <- prediction.epi(theta_ma)
pred <- prediction.epi(theta_la)
# pred <- prediction.epi(new_theta)



  beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM


likelihoods_list <- list()

for (iter in 1:length(beta.params.msm$alpha)) { 
  likelihoods_list[[length(likelihoods_list)+1]] <- 
    data.frame(i = iter, x = seq(0,1,0.01), y = dbeta(x=seq(0,1,0.01), beta.params.msm$alpha[[iter]], beta.params.msm$beta[[iter]]))
}

likelihoods_df <- do.call(rbind.data.frame, likelihoods_list)

ggplot(likelihoods_df, aes(x=x,y=y, group=i)) + 
  geom_line() + 
  geom_vline(data = data.frame(i = 1:length(p.msm.cases), xint = p.msm.cases), mapping = aes(xintercept = xint)) + 
  facet_wrap(~i)




# What are the actual likelihoods being computed here? 


pred <- prediction.epi(theta_ma)
beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM
ma_values <- dbeta(x=p.msm.cases,beta.params.msm$alpha, beta.params.msm$beta, log=TRUE)


pred <- prediction.epi(theta_la)
beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM
la_values <- dbeta(x=p.msm.cases,beta.params.msm$alpha, beta.params.msm$beta, log=TRUE)


pred <- prediction.epi(theta_ma)
beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM
ma_values_div10 <- dbeta(x=p.msm.cases,beta.params.msm$alpha/10, beta.params.msm$beta/10, log=TRUE)


pred <- prediction.epi(theta_la)
beta.params.msm <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.msm"]])), p.msm.cases.var) #beta distn parameter estimates for proportion of males cases who are MSM
la_values_div10 <- dbeta(x=p.msm.cases,beta.params.msm$alpha/10, beta.params.msm$beta/10, log=TRUE)



# What does the p.early target look like 

  beta.params.early <- estBetaParamsVar(as.numeric(unlist(pred[["prev"]][["p.diag.early"]])), p.early.dat$var) #beta distn parameter estimates for proportion of ALL cases (early and late latent) that are early (P,S or early latent)
likelihood_dist <- dbeta(x=seq(0,1,0.01),beta.params.early$alpha, beta.params.early$beta)

ggplot(data.frame(x=seq(0,1,0.01), y=likelihood_dist), aes(x=x,y=y)) + geom_line() + geom_vline(data = data.frame(x = p.early.dat$mean), mapping = aes(xintercept=x))
