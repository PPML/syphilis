

# Setup and Dependencies

devtools::load_all()
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(here)
theme_set(theme_bw())



# Breakdown Function

breakdown_simulation <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

  # rowSum_intervention -- sum up the populations specified by the indices 
	# given (example: diagnosed_females) by year (in the rows). 
	# 
	# Used for Prevalence
	rowSum_intervention <- function(idxs) { 
	  return(rowSums(sol[intervention_period, idxs]))
	}

	data.frame(
		year = intervention_period/52 + model_to_gregorian_difference,
		popsize = rowSums(sol[intervention_period, allpop.index]),


		# Prevalence by Sex
		prevalence = rowSum_intervention(infected.index)
	)
}


#' Simulate Each Intervention and Extract Statistics Using the Passed info_extractor Function
#' 
#' @param theta model parametrization 
#' @param info_extractor a function operating on sol providing a dataframe of summary statistics such 
#'   as prevalent infections, diagnoses, incidence, (by sex, risk, stage, reinfection), etc.
#' 
#' @example 
#' state <- 'LA' # try with 'MA' too!
#' load_globals(model.end = 116) # set the number of years to run out the model to 2050
#' load.start() # parametrize the model
#' theta <- get(paste0('theta_', tolower(state))) # get optimized for the model in the specified state.
#' df <- simulate_interventions(theta, info_extractor = breakdown_simulation)

simulate_interventions <- function(theta, info_extractor) {
	intervention_outcomes <- list()
	for (intervention in c('msm_annual')) {
		e <- constructSimulationEnvironment(theta)
		e$params$output_weekly <- TRUE
		modify_simulation_environment_for_an_intervention(e, intervention)
		runSimulation(e)
		intervention_outcomes[[length(intervention_outcomes)+1]] <- 
			cbind.data.frame(intervention = intervention, info_extractor(e$sim$out))
	}
	do.call(rbind.data.frame, intervention_outcomes)
}


simulate_with_param_replacements <- function(param_list = NULL) { 
	rm(theta_la)
	data(theta_la)

  state <<- 'LA'
	# Simulate and Breakdown Simulation Statistics for Each State and Intervention
	load_globals(model.end = 116) # set the number of years to run out the model to 2050
	load.start() # parametrize the model

	theta_la[param_list] <- theta_ma[param_list]
	theta <- theta_la # get optimized for the model in the specified state.
	la_df <- simulate_interventions(theta, info_extractor = breakdown_simulation) # simulate and summarize each intervention

	return(la_df)
}


base_df <- simulate_with_param_replacements()

partnership_params <- grep("c.min|rp", names(theta_la), value=T)
assortative_params <- grep("theta|epsilon|pi", names(theta_la), value=T)
msm_assortative_params <- grep("theta.4|theta.5|epsilon.4|epsilon.5|pi.msm|theta.hiv", names(theta_la), value=T)
trt_params <- grep("trt", names(theta_la), value=T)

partnership_df <- simulate_with_param_replacements(partnership_params)
assortative_df <- simulate_with_param_replacements(assortative_params)
msm_assortative_df <- simulate_with_param_replacements(msm_assortative_params)
trt_df <- simulate_with_param_replacements(trt_params)

partnership_assortative_df <- simulate_with_param_replacements(c(partnership_params, assortative_params))
assortative_trt_df <- simulate_with_param_replacements(c(assortative_params, trt_params))
trt_partnership_df <- simulate_with_param_replacements(c(trt_params, partnership_params))

all_df <- simulate_with_param_replacements(c(partnership_params, assortative_params, trt_params))

df2 <- rbind.data.frame(
cbind.data.frame(variable = 'base', base_df),
cbind.data.frame(variable = 'partnership', partnership_df),
cbind.data.frame(variable = 'assortative', assortative_df),
cbind.data.frame(variable = 'msm_assortative', assortative_df),
cbind.data.frame(variable = 'trt', trt_df),
cbind.data.frame(variable = 'partnership+assortative', partnership_assortative_df),
cbind.data.frame(variable = 'assortative+trt', assortative_trt_df),
cbind.data.frame(variable = 'trt+partnership', trt_partnership_df),
cbind.data.frame(variable = 'partnership+assortative+trt', all_df))
  
df2 %<>% select(-intervention)

df2 %<>% mutate(prevalence_per_100000 = prevalence / popsize * 100000)

p <- ggplot(df2, aes(x = year, y = prevalence_per_100000, color = variable)) + 
  geom_line()

ggplotly(p)
