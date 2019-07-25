
# Setup and Dependencies

devtools::load_all()
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(tidyr)
library(here)
theme_set(theme_bw())
library('ghibli')


# Breakdown Function

breakdown_simulation <- function(sol) { 

	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)


  # helper functions 

	# rowSum_and_decumulate -- sum up the populations specified by the indices 
	# given (example: diagnosed_females) by year (in the rows) and subtract the 
	# same series lagged by one year (52 months) to ``de-cumulate`` the series. 
	# 
	# Used for Diagnoses and Incidence
	rowSum_and_decumulate <- function(idxs) { 
	  return(rowSums(sol[intervention_period, idxs]) - rowSums(sol[intervention_period - 52, idxs]))
	}

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

		# Incidence by Sex
		incidence = rowSum_and_decumulate(c(inc.index, incr.index)),

		# Prevalence
		prevalence = rowSum_intervention(infected.index),

    # Tests
		tests = rowSum_intervention(tested.index)

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
	for (intervention in c('basecase', 'annual', 'msm_annual',
	'msm_annual_hr_msm_quarterly', 'prior_diagnosis_quarterly',
	'msm_annual_all_hr_annual', 'msm_annual_all_hr_half_annual')) {
		e <- constructSimulationEnvironment(theta)
		e$params$output_weekly <- TRUE
		modify_simulation_environment_for_an_intervention(e, intervention)
		runSimulation(e)
		intervention_outcomes[[length(intervention_outcomes)+1]] <- 
			cbind.data.frame(intervention = intervention, info_extractor(e$sim$out))
	}
	do.call(rbind.data.frame, intervention_outcomes)
}


#### Simulation #####

# Simulate and Breakdown Simulation Statistics for Each State and Intervention
for (state in c('LA', 'MA')) { 
	load_globals(model.end = 116) # set the number of years to run out the model to 2050
	load.start() # parametrize the model
	theta <- get(paste0('theta_', tolower(state))) # get optimized for the model in the specified state.
	assign(paste0(tolower(state), '_df'),
				 simulate_interventions(theta, info_extractor = breakdown_simulation)) # simulate and summarize each intervention
}

df <- rbind.data.frame(
  cbind.data.frame(state = 'Louisiana', la_df),
  cbind.data.frame(state = 'Massachusetts', ma_df))

df %<>% filter(year %in% c(2017, 2027))

df %<>% gather('variable', 'value', -state, -year, -intervention)

# df2 <- reshape2::dcast(setDT(df), state + intervention ~ year, 
#        value.var = c('variable', 'value'))

df %<>% spread('year', 'value')

df$value <- df$`2027` - df$`2017`

df %<>% select(-c('2027', '2017'))

df %<>% spread(variable, value)



basecase_df <- df %>% filter(intervention == 'basecase')

df %<>% filter(intervention != 'basecase')

df %<>% merge(basecase_df, by = c('state'), suffixes = c('', '.basecase'))

df %<>% mutate(additional_tests = tests - tests.basecase, 
               incident_cases_averted = incidence.basecase - incidence,
							 prevalent_cases_averted = prevalence.basecase - prevalence)

df %<>% select(state, intervention, additional_tests, incident_cases_averted, prevalent_cases_averted)

df %<>% mutate(nnt_prevalent = additional_tests / prevalent_cases_averted, nnt_incident = additional_tests / incident_cases_averted) %>% 
  select(state, intervention, nnt_prevalent, nnt_incident)

df2 <- df %>% gather("variable", "value", -state, -intervention)

df2 %<>% mutate(variable = recode(variable, nnt_prevalent = 'prevalent', nnt_incident = 'incident'))

ggplot(df2, aes(x = intervention, y = value, fill = intervention)) + 
  geom_bar(stat = 'identity', alpha = 0.7) + 
	scale_y_log10() + 
	coord_flip() + 
	facet_grid(state~variable) + 
	scale_fill_manual(values = ghibli_palette("MononokeMedium")) + 
	ylab("Number of Tests Needed to Avert One Infection") + 
	xlab("Intervention Scenario") + 
	ggtitle("Number Needed to Test To Avert One Incident/Prevalent Infection", "As Compared to the Base Case in 2027")

ggsave("NumberNeededToTest_Syphilis.png")
