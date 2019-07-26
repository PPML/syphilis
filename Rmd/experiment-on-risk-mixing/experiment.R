# What if Louisiana's Risk Mixing Parameters Were the Same as Massachusetts'

# Let's plot the outcomes (diagnoses, incidence, prevalence)
# by Sex, Stage, and Risk

# Setup and Dependencies

devtools::load_all()
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(here)
theme_set(theme_bw())
library(cowplot)


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

    # Diagnosis by Sex
		diagnosed_women = rowSum_and_decumulate(diagnosed_females),
		diagnosed_msw = rowSum_and_decumulate(diagnosed_msw),
		diagnosed_msm = rowSum_and_decumulate(diagnosed_msm),

    # Diagnosis by Stage
		diagnosed_primary = rowSum_and_decumulate(d1.index),
		diagnosed_secondary = rowSum_and_decumulate(d2.index),
		diagnosed_early_latent = rowSum_and_decumulate(d3.index),
		diagnosed_late_latent = rowSum_and_decumulate(d4.index),

    # Diagnosis by Risk
		diagnosed_high_activity = rowSum_and_decumulate(diagnosed_high_activity),
		diagnosed_low_activity = rowSum_and_decumulate(diagnosed_low_activity),

		# Incidence by Sex
		incidence_women = rowSum_and_decumulate(incidence_females),
		incidence_msm = rowSum_and_decumulate(incidence_msm),
		incidence_msw = rowSum_and_decumulate(incidence_msw),

		# Incidence by Reinfection Status
		incidence_first_time = rowSum_and_decumulate(inc.index),
		incidence_reinfection = rowSum_and_decumulate(incr.index),

		# Incidence by Risk
		incidence_high_activity = rowSum_and_decumulate(incidence_high_activity),
		incidence_low_activity = rowSum_and_decumulate(incidence_low_activity),

		# Prevalence by Sex
		prevalence_women = rowSum_intervention(prevalence_females),
		prevalence_msw = rowSum_intervention(prevalence_msw),
		prevalence_msm = rowSum_intervention(prevalence_msm),

		# Prevalence by Stage
		prevalence_incubating = rowSum_intervention(c(e.index, er.index)),
		prevalence_primary = rowSum_intervention(c(prim.index, primr.index)),
		prevalence_secondary = rowSum_intervention(c(sec.index, secr.index)),
		prevalence_early_latent = rowSum_intervention(c(early.index, earlyr.index)),
		prevalence_late_latent = rowSum_intervention(c(latent.index, latentr.index)),

		# Prevalence by Risk
		prevalence_high_activity = rowSum_intervention(prevalence_high_activity),
		prevalence_low_activity = rowSum_intervention(prevalence_low_activity)
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
	for (intervention in c('basecase', 'annual', 'twelve_times_annual', 'msm_annual',
	'msm_annual_hr_msm_quarterly', 'prior_diagnosis_quarterly',
	'msm_8x_annually', 'twice_annual', 'msm_annual_all_hr_annual'
)) {
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

rm(theta_la)
rm(theta_ma)
data(theta_la)
data(theta_ma)
# risk_mixing_parameters <- grep('epsilon', names(theta_la), value = T)
# theta_la[risk_mixing_parameters] <- theta_ma[risk_mixing_parameters]

# subpopulation_mixing_parameters <- grep('theta', names(theta_la), value = T)
# theta_la[subpopulation_mixing_parameters] <- theta_ma[subpopulation_mixing_parameters]

# treatment_rate_parameters <- grep('trt', names(theta_la), value=T)
# theta_la[treatment_rate_parameters] <- theta_ma[treatment_rate_parameters]

# partnership_parameters <- grep('c.min|rp', names(theta_la), value=T)
# theta_la[partnership_parameters] <- theta_ma[partnership_parameters]

# age_assortative_params <- grep('pi', names(theta_la), value=T)
# theta_la[age_assortative_params] <- theta_ma[age_assortative_params]

# theta_la <- theta_ma

# # Simulate and Breakdown Simulation Statistics for Each State and Intervention
for (state in c('LA', 'MA')) { 
	load_globals(model.end = 116) # set the number of years to run out the model to 2050
	load.start() # parametrize the model
	theta <- get(paste0('theta_', tolower(state))) # get optimized for the model in the specified state.
	assign(paste0(tolower(state), '_df'),
				 simulate_interventions(theta, info_extractor = breakdown_simulation)) # simulate and summarize each intervention
}


#### Data Formatting ####

# Combine Data for Each State
raw_df <- rbind.data.frame(
	cbind.data.frame(state = 'Louisiana', la_df),
	cbind.data.frame(state = 'Massachusetts', ma_df))

# Convert to a Long Format, putting columns of outcomes into the 'variable' column.
# Then separate the variable column into 'outcome' and 'variable'. 
raw_df %>% gather(variable, value, -state, -year, -intervention) %>% 
	separate(col = variable, into = c('outcome', 'variable'), sep = '_', extra = 'merge')-> df

# Pull out the popsize from the basecase for merging back in as a column.
popsize <- filter(df, outcome == 'popsize', intervention == 'basecase') %>% 
  select(-intervention, -outcome, -variable)

# Remove the popsize from the outcomes dataframe.
df %<>% filter(outcome != 'popsize')

# Merge the popsize in as a column.
df %<>% merge(popsize, by = c('state', 'year'), suffixes = c('', 'popsize'))

# Rename the popsize column.
df %<>% rename(popsize = valuepopsize)

# Set outcomes as relative to population size and remove popsize column.
df %<>% mutate(value = value / popsize) %>% select(-popsize)



#### Plotting #### 


outcome_plotter <- function(df, outcome, by, state) { 

  stopifnot(outcome %in% c('prevalence', 'diagnosed', 'incidence'))
  stopifnot(by %in% c('sex', 'stage', 'risk', 'reinfection'))
	stopifnot(state %in% c('Louisiana', 'Massachusetts'))

	state <- enquo(state)
	outcome_str <- outcome
	outcome <- enquo(outcome)

  var_set <- switch(by,
	  sex = c('women', 'msm', 'msw'),
		stage = c('incubating', 'primary', 'secondary', 'early_latent', 'late_latent'),
		risk = c('high_activity', 'low_activity'),
		reinfection = c('first_time', 'reinfection'))

  formatted_outcome <- switch(outcome_str,
	  prevalence = 'Prevalent Infections ',
		incidence = 'Incident Infections ',
		diagnosed = 'Diagnosed Infections ')

	title <- paste0(
	  formatted_outcome,
		"per 100,000 ",
		if (outcome_str == 'prevalence') '' else 'per Year ',
		'by ',
		tools::toTitleCase(by))

	color_palette <- switch(outcome_str,
	  prevalence = 'Reds',
		diagnosed = 'Greens',
		incidence = 'Blues')

  ggplot(filter(df, state == !!state, outcome == !!outcome, variable %in% var_set), aes(x = year, y = value * 100000, fill = variable)) + 
	geom_area(colour="black", size=.2, alpha=.4) + 
	scale_fill_brewer(palette=color_palette) + 
	xlab("Year") + 
	facet_wrap(~ intervention) + 
	ylab("Prevalent Infections per 100,000") + 
	scale_x_continuous(breaks = seq(2016, 2027, by = 2)) + 
	theme(axis.text.x = element_text(angle = 90)) + 
	ggtitle(title, subtitle = state)

}


combinations_table <- expand.grid(
  state = c('Louisiana', 'Massachusetts'),
	outcome = c('prevalence', 'diagnosed', 'incidence'),
	by = c('sex', 'stage', 'risk', 'reinfection'))


combinations_table %<>% 
  filter(! (outcome == 'incidence' & by == 'stage'))

combinations_table %<>% 
  filter(! (outcome %in% c('prevalence', 'diagnosed') & by == 'reinfection'))

combinations_table %<>% 
  mutate_if(is.factor, as.character)

shortnames <- setNames(nm = statenames, object = names(statenames))

for (row_i in 1:nrow(combinations_table)) { 

  state <- combinations_table[row_i, 'state']
  outcome <- combinations_table[row_i, 'outcome']
  by <- combinations_table[row_i, 'by']

	state_abb <- tolower(shortnames[[state]])

	plot_name <- 
	  paste0(
		outcome, '_', 'by_', by, '_', state_abb
		)

  plt <- outcome_plotter(df, outcome = outcome, by = by, state = state)

  assign(plot_name, plt)

	# ggsave(filename = paste0('individual_plots/', plot_name, '.png'), plot = plt)
}


# Combine Plots
for (state in c('la', 'ma')) { 
	plt <- cowplot::plot_grid(
		get(paste0('prevalence_by_sex_', state)),
		get(paste0('incidence_by_sex_', state)),
		get(paste0('diagnosed_by_sex_', state)),
		get(paste0('prevalence_by_stage_', state)),
		get(paste0('incidence_by_reinfection_', state)),
		get(paste0('diagnosed_by_stage_', state)),
		get(paste0('prevalence_by_risk_', state)),
		get(paste0('incidence_by_risk_', state)),
		get(paste0('diagnosed_by_risk_', state)),
		ncol = 3)
	ggsave(plot = plt, paste0('outcomes_', state, '_overwrite_risk_mixing', '.svg'), height = 30, width = 30, dpi = 600)
}
