#' Modify the Screening Matrix to Implement Targeted Testing and Treatment Interventions
#' 
#' Check that the intervention is defined appropriately (either not at all or
#' as basecase, or as one of the interventions in the interventions$codename
#' column). Then use a switch on the intervention to call the
#' adjust_screening_for_intervention function with the right arguments.
#' 
modify_simulation_environment_for_an_intervention <- function(e, intervention) {

	if (intervention == 'basecase') return(e)

	# if the intervention is not among those defined in interventions, error
	if (! intervention %in% interventions$codename) {
		stop("the intervention must be defined in interventions definition dataframe")
	}

  # turn the intervention codename into arguments for adjust_screening_for_intervention
	switch(intervention, 
									annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 1),
						twice_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 2),
							msm_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 1),
				msm_twice_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 2),
			 msm_hivpos_annual = adjust_screening_for_intervention(e, pop_indices = m5, new_level = 1),
 msm_hivpos_twice_annual = adjust_screening_for_intervention(e, pop_indices = m5, new_level = 2),
			 msm_hivneg_annual = adjust_screening_for_intervention(e, pop_indices = m4, new_level = 1),
 msm_hivneg_twice_annual = adjust_screening_for_intervention(e, pop_indices = m4, new_level = 2),
	)

  return(invisible(NULL))
}

#' Modify Screening Matrix for A Targeted Screening Intervention
#' 
#' Using the intervention_years as rows and the pop_indices as columns, 
#' select those entries in the screening matrix during the intervention 
#' time period and among the targeted population which are less than 
#' the new_level of screening and update them to the new_level.
#'
adjust_screening_for_intervention <- function(e, intervention_start = 105, pop_indices, new_level) {

	e$params[['output_weekly']] <- TRUE

	# First check that the intervention years are within the modeled time period
	# and that the time varying parameters extend to at least the end of the
	# intervention period.
	max_year <- nrow(e$screen)
	intervention_years <- intervention_start:max_year

	if (any(c(e$model.end, nrow(e$screen), length(e$behav),
							nrow(e$rep.trend)) < max_year )) {
		stop("The model.end needs to be at or after the end of the intervention period.
Additionally, it must be defined before constructing the simulation environment
so that the screening, transmission probability, and reporting probability
time trends are constructed with entries through at least the end of the 
intervention period.")
	}

	# Determine which entries during the intervention_years among the 
	# selected populations (pop_indices) are below the new_level of 
	# screening and then replace those entries' values with the 
	# new_level.
	intervention_timesteps = seq(from=min(intervention_years)*52, to=max(intervention_years))
	need_increase <- e$screen[intervention_years, pop_indices] < new_level
	e$screen[intervention_years, pop_indices][need_increase] <- new_level
  
	# Remember that the version of screen the C++ uses is the one from 
	# the params list, so be sure to update that as well.
	e$params$screen <- e$screen

	return(invisible(NULL))
}


#' Extract Data from Simulation Outcomes for Screening Intervention Analysis
#' 
#' In our screening intervention analysis we will calculate the following 
#' values: 
#' 
#' - Number of Additional Tests
#' - Number of Incident Cases 
#' - Number of Prevalent Cases
#' - Number Needed to Test to Avert One Additional Incident Case
#' - Number Needed to Test to Avert One Prevalent Case
#' 
#' In order to calculate these, we will extract the following data from a
#' simulation: prevalent_infections, incidents, tests (per year).
#' 
#' @example 
#' devtools::load_all()
#' state <- c('LA', 'MA')[sample.int(2,1)]
#' load.start()
#' theta <- get(paste0('theta_', tolower(state)))
#' e <- constructSimulationEnvironment(theta)
#' intervention <- interventions$codename[sample.int(length(interventions$codename),1)]
#' modify_simulation_environment_for_an_intervention(e, intervention = intervention)
#' runSimulation(e)
#' df <- extractInterventionStatistics(e, intervention = intervention)
extractInterventionStatistics <- function(e, intervention) {
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

	data.frame(
		intervention = intervention,
		year = intervention_period/52 + model_to_gregorian_difference,
		prevalent_infections = rowSums(e$sim$out[1+intervention_period, 1+infected.index]),
		incidents_per_year = rowSums(e$sim$out[1+intervention_period, 1+c(inc.index, incr.index)]) - 
			rowSums(e$sim$out[intervention_period, 1+c(inc.index, incr.index)]),
		tests_per_year = rowSums(e$sim$out[1+intervention_period,1+tested.index]) - 
			rowSums(e$sim$out[intervention_period,1+tested.index]) 
		)
}


#' Simulate Each Intervention with a Given Parametrization
#' 
#' Go through each of the interventions (including the basecase),
#' simulate them, extract the Intervention Statistics, and 
#' construct a dataframe with the results.
#'
#' @example 
#' devtools::load_all()
#' state <- c('LA', 'MA')[sample.int(2,1)]
#' load.start()
#' theta <- get(paste0('theta_', tolower(state)))
#' df <- simulate_interventions(theta)
simulate_interventions <- function(theta) {
	intervention_outcomes <- list()
	for (intervention in c('basecase', interventions$codename)) {
		e <- constructSimulationEnvironment(theta)
		modify_simulation_environment_for_an_intervention(e, intervention)
		runSimulation(e)
		intervention_outcomes[[length(intervention_outcomes)+1]] <- 
			extractInterventionStatistics(e, intervention)
	}
	do.call(rbind.data.frame, intervention_outcomes)
}


#' Format the Intervention Statistics 
#' 
#' Take the output from simulate_interventions do the following operations: 
#' - Calculate Additional Tests
#' - Calculate Change In Prevalent Infections 
#' - Calculate Change in 
#' 
#' @param df A dataframe resulting from simulate_interventions.
#'
#' @example 
#' devtools::load_all()
#' state <- c('LA', 'MA')[sample.int(2,1)]
#' load.start()
#' theta <- get(paste0('theta_', tolower(state)))
#' df <- simulate_interventions(theta)
#' df <- format_intervention_statistics(df)
#' 
#' # prevalent infections in each
#' ggplot(df, aes(x = year, y = prevalent_infections, color = intervention)) + 
#'   geom_line()
#'
#' # prevalent infections averted
#' ggplot(df, aes(x = year, y = -1 * chg_prevalent_infs, color = intervention)) + 
#'   geom_line() +
#'   scale_y_continuous(trans = symlog_trans())
#' 
#' # prevalent infections averted in 2021
#' ggplot(filter(df, year == 2021), aes(x=intervention, y = -chg_prevalent_infs, fill = intervention, color = intervention)) + 
#'	 geom_bar(stat='identity', alpha = 0.6) + 
#'	 geom_text(aes(x = intervention, label = intervention, y = -chg_prevalent_infs*.1 + .5), color = 'black', alpha = 0.5, size=4) + # , angle = 25) + 
#'	 scale_y_continuous(trans = symlog_trans()) + 
#'	 coord_flip() + 
#'	 theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = 'none')
#'
#' 
#' # incident infections in each
#' ggplot(df, aes(x = year, y = incidents, color = intervention)) + 
#'   geom_line()
#'
#' # prevalent infections averted
#' ggplot(df, aes(x = year, y = -1 * chg_prevalent_infs, color = intervention)) + 
#'   geom_line() +
#'   scale_y_continuous(trans = symlog_trans())
#' 
#' # prevalent infections averted in 2021
#' ggplot(filter(df, year == 2021), aes(x=intervention, y = -chg_prevalent_infs, fill = intervention, color = intervention)) + 
#'	 geom_bar(stat='identity', alpha = 0.6) + 
#'	 geom_text(aes(x = intervention, label = intervention, y = -chg_prevalent_infs*.1 + .5), color = 'black', alpha = 0.5, size=4) + # , angle = 25) + 
#'	 scale_y_continuous(trans = symlog_trans()) + 
#'	 coord_flip() + 
#'	 theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = 'none')
format_intervention_statistics <- function(df) {
	basecase <- df %>% dplyr::filter(intervention == 'basecase') %>% dplyr::select(-intervention)
	df <- merge(df, basecase, by=c('year'), suffixes = c('', '.basecase'))
	df <- df %>% dplyr::mutate(
	    chg_prevalent_infs = prevalent_infections - prevalent_infections.basecase,
			chg_prevalent_infs_pct = chg_prevalent_infs / prevalent_infections.basecase,
			chg_incident_infs = incidents - incidents.basecase,
			chg_incident_infs_pct = chg_incident_infs / incidents.basecase,
			chg_tests = tests_per_year - tests_per_year.basecase,
			chg_tests_pct = chg_tests / tests_per_year.basecase,
			nnt_prevalent = -chg_prevalent_infs / chg_tests,
			nnt_incident = -chg_incident_infs / chg_tests
	  )
	df <- dplyr::select(df, -c('prevalent_infections.basecase', 'incidents.basecase', 'tests_per_year.basecase'))
	df <- dplyr::arrange(df, intervention, year)
	return(df)
}
