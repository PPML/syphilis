#' Modify the Screening Matrix to Implement Targeted Testing and Treatment Interventions
#' 
#' Check that the intervention is defined appropriately (either not at all or
#' as basecase, or as one of the interventions in the interventions$codename
#' column). Then use a switch on the intervention to call the
#' adjust_screening_for_intervention function with the right arguments.
#' 
modify_simulation_environment_for_an_intervention <- function(e) {

  # return e unmodified if no intervention is specified
  if (! exists('intervention') || (exists('intervention') && intervention == 'basecase')) return(e)

  # if the intervention is not among those defined in interventions, error
	if (exists('intervention') && 
	    ! intervention %in% interventions$codename) {
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
		msm_hiv_twice_annual = adjust_screening_for_intervention(e, pop_indices = m4, new_level = 2),
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
adjust_screening_for_intervention <- function(e, intervention_years = 105:109, pop_indices, new_level) {

	# First check that the intervention years are within the modeled time period
	# and that the time varying parameters extend to at least the end of the
	# intervention period.
	max_year <- max(intervention_years)

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
extractInterventionStatistics <- function(e, intervention) {
	data.frame(
		intervention = intervention,
		year = 2017:2021,
		prevalent_infections = rowSums(e$sim$out[1+intervention_years, 1+infected.index]),
		incidents = rowSums(e$sim$out[1+intervention_years, 1+c(inc.index, incr.index)]) - 
			rowSums(e$sim$out[intervention_years, 1+c(inc.index, incr.index)]),
		tests = rowSums(e$sim$out[1+intervention_years,1+tested.index]) - 
			rowSums(e$sim$out[intervention_years,1+tested.index]) 
		)
}


#' Simulate Each Intervention with a Given Parametrization
#' 
#' Go through each of the interventions (including the basecase),
#' simulate them, extract the Intervention Statistics, and 
#' construct a dataframe with the results.
simulate_interventions <- function(theta) {
	intervention_outcomes <- list()
	for (intervention in c('basecase', interventions$codename)) {
		e <- constructSimulationEnvironment(theta)
		modify_simulation_environment_for_an_intervention(e)
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
format_intervention_statistics <- function(df) {
}
