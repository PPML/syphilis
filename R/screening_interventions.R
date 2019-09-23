#' Modify the Screening Matrix to Implement Targeted Testing and Treatment Interventions
#' 
#' Check that the intervention is defined appropriately (either not at all or
#' as basecase, or as one of the interventions in the interventions$codename
#' column). Then use a switch on the intervention to call the
#' adjust_screening_for_intervention function with the right arguments.
#' 
modify_simulation_environment_for_an_intervention <- function(e, intervention) {

	if (intervention == 'basecase') return(e)

  first_intervention_year <- min(intervention_years) 

  e$params$screen[(first_intervention_year*52):nrow(e$params$screen),setdiff(1:40, y.f)] <- 0
  e$params$screen_repeat[(first_intervention_year*52):nrow(e$params$screen_repeat),setdiff(1:40, y.f)] <- 0

  # turn the intervention codename into arguments for adjust_screening_for_intervention
	switch(intervention, 
									annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 1),
						twice_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 2),
			 twelve_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 12),
   twenty_four_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 24),
   sixteen_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 16),
							msm_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 1),
				msm_twice_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 2),
				msm_quarterly = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 4),
			 msm_hivpos_annual = adjust_screening_for_intervention(e, pop_indices = m5, new_level = 1),
   msm_hivpos_twice_annual = adjust_screening_for_intervention(e, pop_indices = m5, new_level = 2),
			 msm_hivneg_annual = adjust_screening_for_intervention(e, pop_indices = m4, new_level = 1),
   msm_hivneg_twice_annual = adjust_screening_for_intervention(e, pop_indices = m4, new_level = 2),
   msm_annual_hr_msm_quarterly = {
		adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 1)
		adjust_screening_for_intervention(
			e, 
			pop_indices = # high sexual activity msm
				setdiff(c(m4, m5), setdiff(1:40, pop$index[pop$k == 'high'])),
			new_level = 4)
	 },
   msm_8x_annually = adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 8),
   hr_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 1),
   hr_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 1),
   msm_annual_hr_and_prior_diagnosis_quarterly = {
     adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 1)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4)
      e$params$screen_repeat[first_intervention_year:nrow(e$params$screen_repeat), ] <- 4
     },

   msm_hr_and_prior_diagnosis_quarterly = {
     adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 4)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4)
      e$params$screen_repeat[first_intervention_year:nrow(e$params$screen_repeat), ] <- 4
     },
   msm_annual_all_hr_annual = {
     adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 1)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 1)
     },
   msm_annual_all_hr_half_annual = {
     adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 1)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = .5)
     },
   prior_diagnosis_quarterly = {
      e$params$screen_repeat[(first_intervention_year*52):nrow(e$params$screen_repeat), ] <- 
        4
     },
   prior_diagnosis_annual = {
      e$params$screen_repeat[(first_intervention_year*52):nrow(e$params$screen_repeat), ] <- 
        1
     },
   high_activity_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 1),
   high_activity_twice_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 2),
   high_activity_quarterly = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4),
   annual_hr_quarterly = { 
     adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 1)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4)
   },
   half_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 0.5),
   quarterly = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 4),
   eight_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 8)
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
adjust_screening_for_intervention <- function(e, intervention_start, pop_indices, new_level, repeat_level = NA) {

  if (missing(intervention_start)) intervention_start <- min(intervention_years) # +1 for C++ indexing of screen

	# First check that the intervention years are within the modeled time period
	# and that the time varying parameters extend to at least the end of the
	# intervention period.
	max_year <- nrow(e$screen)
	intervention_years <- (intervention_start*52):max_year

  
	# if (any(c(e$model.end, nrow(e$screen), # length(e$behav),
	# 						nrow(e$rep.trend)) < max_year )) {
	# 	stop("The model.end needs to be at or after the end of the intervention period.
# Additionally, it must be defined before constructing the simulation environment
# so that the screening, transmission probability, and reporting probability
# time trends are constructed with entries through at least the end of the 
# intervention period.")
	# }

	# Determine which entries during the intervention_years among the 
	# selected populations (pop_indices) are below the new_level of 
	# screening and then replace those entries' values with the 
	# new_level.
	intervention_timesteps = seq(from=min(intervention_years), to=max(intervention_years))
	# need_increase <- e$screen[intervention_years, pop_indices] < new_level
	e$screen[intervention_years, pop_indices] <- new_level
  
	# Remember that the version of screen the C++ uses is the one from 
	# the params list, so be sure to update that as well.
	e$params$screen <- e$screen
	e$params$screen_repeat <- e$params$screen

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
# extractInterventionStatistics <- function(sol) {
# 	intervention_start <- min(intervention_years)
# 	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

# 	data.frame(
# 		year = intervention_period/52 + model_to_gregorian_difference,
# 		prevalent_infections = rowSums(sol[1+intervention_period, 1+infected.index]),
# 		incidents_per_year = rowSums(sol[1+intervention_period, 1+c(inc.index, incr.index)]) - 
# 			rowSums(sol[1 + intervention_period - 52, 1+c(inc.index, incr.index)]),
# 		tests_per_year = rowSums(sol[1+intervention_period,1+tested.index]) - 
# 			rowSums(sol[intervention_period,1+tested.index]),
# 		popsize = rowSums(sol[1+intervention_period, 1+allpop.index])
# 		)
# }


extract_outcomes <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, (intervention_start + cal.period)*52)

  decumulate <- function(idxs) { 
    return(
      rowSums(sol[intervention_period, 1 + idxs]) - 
      rowSums(sol[intervention_period - 52, 1 + idxs])
    )} 

  retrieve <- function(idxs) { 
    rowSums(sol[intervention_period, 1 + idxs])
  }

  data.frame(
    year = intervention_period/52 + model_to_gregorian_difference,

    # popsizes
    popsize_all = retrieve(allpop.index),
    # popsize_young = retrieve(all_young),
    # popsize_old = retrieve(all_old),

    # popsize_females = retrieve(all_females),
    popsize_young_females = retrieve(all_young_females),
    popsize_old_females = retrieve(all_old_females),

    # popsize_msw = retrieve(all_msw),
    # popsize_young_msw = retrieve(all_young_msw),
    # popsize_old_msw = retrieve(all_old_msw),

    # popsize_msm = retrieve(all_msm),
    # popsize_young_msm = retrieve(all_young_msm),
    # popsize_old_msm = retrieve(all_old_msm),

    # popsize_males = retrieve(all_males),
    popsize_young_males = retrieve(all_young_males),
    popsize_old_males = retrieve(all_old_males),

    # diagnoses 
    diagnosed_all = decumulate(c(d1.index, d2.index, d3.index)), # diagnosed early syphilis
    # diagnosed_young = decumulate(diagnosed_young), 
    # diagnosed_old = decumulate(diagnosed_old), 

    # diagnosed_females = decumulate(diagnosed_females),
    diagnosed_young_females = decumulate(diagnosed_young_females),
    diagnosed_old_females = decumulate(diagnosed_old_females),

    # diagnosed_msw = decumulate(diagnosed_msw),
    # diagnosed_young_msw = decumulate(diagnosed_young_msw),
    # diagnosed_old_msw = decumulate(diagnosed_old_msw),

    # diagnosed_msm = decumulate(diagnosed_msm),
    # diagnosed_young_msm = decumulate(diagnosed_young_msm),
    # diagnosed_old_msm = decumulate(diagnosed_old_msm),

    # diagnosed_males = decumulate(diagnosed_males),
    diagnosed_young_males = decumulate(diagnosed_young_males),
    diagnosed_old_males = decumulate(diagnosed_old_males),

    # diagnosed_reinfected = decumulate(dr.index),

    # incidence
    incidence_all = decumulate(c(inc.index, incr.index)),
    # incidence_young = decumulate(incidence_young),
    # incidence_old = decumulate(incidence_old),

    # incidence_females = decumulate(incidence_females),
    incidence_young_females = decumulate(incidence_young_females),
    incidence_old_females = decumulate(incidence_old_females),

    # incidence_msw = decumulate(incidence_msw),
    # incidence_young_msw = decumulate(incidence_young_msw),
    # incidence_old_msw = decumulate(incidence_old_msw),

    # incidence_msm = decumulate(incidence_msm),
    # incidence_young_msm = decumulate(incidence_young_msm),
    # incidence_old_msm = decumulate(incidence_old_msm),

    # incidence_males = decumulate(incidence_males),
    incidence_young_males = decumulate(incidence_young_males),
    incidence_old_males = decumulate(incidence_old_males),

    # incidence_reinfected = decumulate(incr.index),
    # incidence_new = decumulate(inc.index),
    # incidence_high_activity = decumulate(incidence_high_activity),
    # incidence_low_activity = decumulate(incidence_low_activity),



    # prevalence
    prevalence_all = retrieve(early_infected_index),
    # prevalence_young = retrieve(prevalence_early_young),
    # prevalence_old = retrieve(prevalence_early_old),

    # prevalence_females = retrieve(prevalence_early_females),
    prevalence_young_females = retrieve(prevalence_early_young_females),
    prevalence_old_females = retrieve(prevalence_early_old_females),

    # prevalence_msw = retrieve(prevalence_early_msw),
    # prevalence_young_msw = retrieve(prevalence_early_young_msw),
    # prevalence_old_msw = retrieve(prevalence_early_old_msw),

    # prevalence_msm = retrieve(prevalence_early_msm),
    # prevalence_young_msm = retrieve(prevalence_early_young_msm),
    # prevalence_old_msm = retrieve(prevalence_early_old_msm),

    # prevalence_males = retrieve(prevalence_early_males),
    prevalence_young_males = retrieve(prevalence_early_young_males),
    prevalence_old_males = retrieve(prevalence_early_old_males)

    # prevalence_reinfected = retrieve(reinfected.index),
    # prevalence_new = retrieve(newly_infected.index),
    # prevalence_high_activity = retrieve(prevalence_high_activity),
    # prevalence_low_activity = retrieve(prevalence_low_activity),
    # prevalence_noninfectious = retrieve(noninfectious_index),
    # prevalence_infectious = retrieve(infectious_index),

    # tests
    # tests_all = decumulate(tested.index)
  )

}


# extract_diagnoses_by_sex <- function(sol) { 
# 	intervention_start <- min(intervention_years)
# 	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

#   decumulate <- function(idxs) { 
#     return(
#       rowSums(sol[intervention_period, idxs]) - 
#       rowSums(sol[intervention_period - 52, idxs])
#     )} 

#   data.frame(
#     year = intervention_period/52 + model_to_gregorian_difference,
#     diagnosed_females = decumulate(diagnosed_females),
#     diagnosed_msw = decumulate(diagnosed_msw),
#     diagnosed_msm = decumulate(diagnosed_msm),
#     popsize_females = rowSums(sol[intervention_period, females]),
#     popsize_msw = rowSums(sol[intervention_period, msw]),
#     popsize_msm = rowSums(sol[intervention_period, msm])
#   )
# }


# extract_incidence_by_sex <- function(sol) { 
# 	intervention_start <- min(intervention_years)
# 	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

#   decumulate <- function(idxs) { 
#     return(
#       rowSums(sol[intervention_period, idxs]) - 
#       rowSums(sol[intervention_period - 52, idxs])
#     )} 

#   data.frame(
#     year = intervention_period/52 + model_to_gregorian_difference,
#     incidence_females = decumulate(incidence_females),
#     incidence_msw = decumulate(incidence_msw),
#     incidence_msm = decumulate(incidence_msm),
#     popsize_females = rowSums(sol[intervention_period, females]),
#     popsize_msw = rowSums(sol[intervention_period, msw]),
#     popsize_msm = rowSums(sol[intervention_period, msm])
#   )
# }


# extract_prevalence_by_sex <- function(sol) { 
# 	intervention_start <- min(intervention_years)
# 	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

#   decumulate <- function(idxs) { 
#     return(
#       rowSums(sol[intervention_period, idxs]) - 
#       rowSums(sol[intervention_period - 52, idxs])
#     )} 

#   data.frame(
#     year = intervention_period/52 + model_to_gregorian_difference,
#     prevalence_females = rowSums(sol[intervention_period, prevalence_females]),
#     prevalence_msw = rowSums(sol[intervention_period, prevalence_msw]),
#     prevalence_msm = rowSums(sol[intervention_period, prevalence_msm]),
#     popsize_females = rowSums(sol[intervention_period, females]),
#     popsize_msw = rowSums(sol[intervention_period, msw]),
#     popsize_msm = rowSums(sol[intervention_period, msm])
#   )
# }


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
simulate_interventions <- function(
  theta, 
  interventions_list = c(
  'msm_annual_hr_msm_quarterly',
  'msm_quarterly',
  'prior_diagnosis_quarterly', 
  'high_activity_quarterly'
  ),
  info_extractor
) {
	intervention_outcomes <- list()
	for (intervention in c('basecase', interventions_list)) {
		e <- constructSimulationEnvironment(theta)
		modify_simulation_environment_for_an_intervention(e, intervention)
		e$params$output_weekly <- TRUE
		runSimulation(e)
		intervention_outcomes[[length(intervention_outcomes)+1]] <- 
			cbind.data.frame(intervention = intervention, info_extractor(e$sim$out))
	}
	do.call(rbind.data.frame, intervention_outcomes)
}

simulate_interventions_trace <- function(
  trace, 
  interventions_list = c(
  'msm_annual_hr_msm_quarterly',
  'msm_quarterly',
  'prior_diagnosis_quarterly', 
  'high_activity_quarterly'
  ),
  info_extractor,
  n
) {
  outcomes <- list()
  if (missing(n)) n <- nrow(trace) # default n to nrow of trace

  for (iter in 1:n) { 
    theta <- unlist(trace[iter,]) # get theta from trace
		intervention_outcomes <- list() # setup list for storing intervention outcomes for theta
		for (intervention in c('basecase', interventions_list)) { # simulate each intervention
			e <- constructSimulationEnvironment(theta)
			modify_simulation_environment_for_an_intervention(e, intervention)
			e$params$output_weekly <- TRUE
			runSimulation(e)
			intervention_outcomes[[length(intervention_outcomes)+1]] <- 
				cbind.data.frame(intervention = intervention, info_extractor(e$sim$out))
		}
		outcomes[[iter]] <- cbind.data.frame(sim = iter, do.call(rbind.data.frame, intervention_outcomes)) # store the trace[iter,] intervention outcomes
	}
  outcomes <- do.call(rbind.data.frame, outcomes) # convert outcomes from list of dfs to df
  return(outcomes)
}


compute_confidence_intervals_for_trace_interventions <- function() { 
  # read in raw data
  la_df <- readRDS(system.file('interventions/la_interventions_df.rds', package='syphilis'))
  ma_df <- readRDS(system.file('interventions/ma_interventions_df.rds', package='syphilis'))

  # > colnames(la_df)
  #  [1] "sim"                      "intervention"
  #  [3] "year"                     "popsize_young_females"
  #  [5] "popsize_old_females"      "popsize_young_males"
  #  [7] "popsize_old_males"        "diagnosed_young_females"
  #  [9] "diagnosed_old_females"    "diagnosed_young_males"
  # [11] "diagnosed_old_males"      "incidence_young_females"
  # [13] "incidence_old_females"    "incidence_young_males"
  # [15] "incidence_old_males"      "prevalence_young_females"
  # [17] "prevalence_old_females"   "prevalence_young_males"
  # [19] "prevalence_old_males"

  # combine the data frames
  la_df <- cbind.data.frame(state = 'Louisiana', la_df)
  ma_df <- cbind.data.frame(state = 'Massachusetts', ma_df)
  df <- rbind.data.frame(la_df, ma_df)

  # group by state, year, intervention, and produce confidence intervals
  df %>% group_by(state, year, intervention) %>% 
    summarize(
      # preserve population sizes
      popsize_young_females = median(popsize_young_females, na.rm=T),
      popsize_old_females = median(popsize_old_females, na.rm=T),
      popsize_young_males = median(popsize_young_males, na.rm=T),
      popsize_old_males = median(popsize_old_males, na.rm=T),
      
      # yearly diagnosed counts -- median, ci_high, and ci_low (95% confidence intervals)
      diagnosed_young_females.median = median(diagnosed_young_females, na.rm=T),
      diagnosed_young_females.ci_high = quantile(diagnosed_young_females, .975, na.rm=T),
      diagnosed_young_females.ci_low = quantile(diagnosed_young_females, .025, na.rm=T),

      diagnosed_old_females.median = median(diagnosed_old_females, na.rm=T),
      diagnosed_old_females.ci_high = quantile(diagnosed_old_females, .975, na.rm=T),
      diagnosed_old_females.ci_low = quantile(diagnosed_old_females, .025, na.rm=T),

      diagnosed_young_males.median = median(diagnosed_young_males, na.rm=T),
      diagnosed_young_males.ci_high = quantile(diagnosed_young_males, .975, na.rm=T),
      diagnosed_young_males.ci_low = quantile(diagnosed_young_males, .025, na.rm=T),

      diagnosed_old_males.median = median(diagnosed_old_males, na.rm=T),
      diagnosed_old_males.ci_high = quantile(diagnosed_old_males, .975, na.rm=T),
      diagnosed_old_males.ci_low = quantile(diagnosed_old_males, .025, na.rm=T),

      # yearly incidence counts -- median, ci_high, and ci_low (95% confidence intervals)
      incidence_young_females.median = median(incidence_young_females, na.rm=T),
      incidence_young_females.ci_high = quantile(incidence_young_females, .975, na.rm=T),
      incidence_young_females.ci_low = quantile(incidence_young_females, .025, na.rm=T),

      incidence_old_females.median = median(incidence_old_females, na.rm=T),
      incidence_old_females.ci_high = quantile(incidence_old_females, .975, na.rm=T),
      incidence_old_females.ci_low = quantile(incidence_old_females, .025, na.rm=T),

      incidence_young_males.median = median(incidence_young_males, na.rm=T),
      incidence_young_males.ci_high = quantile(incidence_young_males, .975, na.rm=T),
      incidence_young_males.ci_low = quantile(incidence_young_males, .025, na.rm=T),

      incidence_old_males.median = median(incidence_old_males, na.rm=T),
      incidence_old_males.ci_high = quantile(incidence_old_males, .975, na.rm=T),
      incidence_old_males.ci_low = quantile(incidence_old_males, .025, na.rm=T),

      # prevalence counts -- median, ci_high, and ci_low (95% confidence intervals)
      prevalence_young_females.median = median(prevalence_young_females, na.rm=T),
      prevalence_young_females.ci_high = quantile(prevalence_young_females, .975, na.rm=T),
      prevalence_young_females.ci_low = quantile(prevalence_young_females, .025, na.rm=T),

      prevalence_old_females.median = median(prevalence_old_females, na.rm=T),
      prevalence_old_females.ci_high = quantile(prevalence_old_females, .975, na.rm=T),
      prevalence_old_females.ci_low = quantile(prevalence_old_females, .025, na.rm=T),

      prevalence_young_males.median = median(prevalence_young_males, na.rm=T),
      prevalence_young_males.ci_high = quantile(prevalence_young_males, .975, na.rm=T),
      prevalence_young_males.ci_low = quantile(prevalence_young_males, .025, na.rm=T),

      prevalence_old_males.median = median(prevalence_old_males, na.rm=T),
      prevalence_old_males.ci_high = quantile(prevalence_old_males, .975, na.rm=T),
      prevalence_old_males.ci_low = quantile(prevalence_old_males, .025, na.rm=T)
    ) -> df2

    saveRDS(df2, file.path(system.file('interventions/', package='syphilis'), 'intvs_summarized.rds'))

  df2 %>% 
  # reshape data into rates per 100k
  dplyr::mutate(
    # prevalence 
    prevalence_youngfemales.ci_high = prevalence_young_females.ci_high / popsize_young_females * 1e5,
    prevalence_oldfemales.ci_high = prevalence_old_females.ci_high / popsize_old_females * 1e5,
    prevalence_youngmales.ci_high = prevalence_young_males.ci_high / popsize_young_males * 1e5,
    prevalence_oldmales.ci_high = prevalence_old_males.ci_high / popsize_old_males * 1e5,

    prevalence_youngfemales.ci_low = prevalence_young_females.ci_low / popsize_young_females * 1e5,
    prevalence_oldfemales.ci_low = prevalence_old_females.ci_low / popsize_old_females * 1e5,
    prevalence_youngmales.ci_low = prevalence_young_males.ci_low / popsize_young_males * 1e5,
    prevalence_oldmales.ci_low = prevalence_old_males.ci_low / popsize_old_males * 1e5,

    prevalence_youngfemales.median = prevalence_young_females.median / popsize_young_females * 1e5,
    prevalence_oldfemales.median = prevalence_old_females.median / popsize_old_females * 1e5,
    prevalence_youngmales.median = prevalence_young_males.median / popsize_young_males * 1e5,
    prevalence_oldmales.median = prevalence_old_males.median / popsize_old_males * 1e5,

    # incidence
    incidence_youngfemales.ci_high = incidence_young_females.ci_high / popsize_young_females * 1e5,
    incidence_oldfemales.ci_high = incidence_old_females.ci_high / popsize_old_females * 1e5,
    incidence_youngmales.ci_high = incidence_young_males.ci_high / popsize_young_males * 1e5,
    incidence_oldmales.ci_high = incidence_old_males.ci_high / popsize_old_males * 1e5,

    incidence_youngfemales.ci_low = incidence_young_females.ci_low / popsize_young_females * 1e5,
    incidence_oldfemales.ci_low = incidence_old_females.ci_low / popsize_old_females * 1e5,
    incidence_youngmales.ci_low = incidence_young_males.ci_low / popsize_young_males * 1e5,
    incidence_oldmales.ci_low = incidence_old_males.ci_low / popsize_old_males * 1e5,

    incidence_youngfemales.median = incidence_young_females.median / popsize_young_females * 1e5,
    incidence_oldfemales.median = incidence_old_females.median / popsize_old_females * 1e5,
    incidence_youngmales.median = incidence_young_males.median / popsize_young_males * 1e5,
    incidence_oldmales.median = incidence_old_males.median / popsize_old_males * 1e5,

    # diagnoses
    diagnosed_youngfemales.ci_high = diagnosed_young_females.ci_high / popsize_young_females * 1e5,
    diagnosed_oldfemales.ci_high = diagnosed_old_females.ci_high / popsize_old_females * 1e5,
    diagnosed_youngmales.ci_high = diagnosed_young_males.ci_high / popsize_young_males * 1e5,
    diagnosed_oldmales.ci_high = diagnosed_old_males.ci_high / popsize_old_males * 1e5,

    diagnosed_youngfemales.ci_low = diagnosed_young_females.ci_low / popsize_young_females * 1e5,
    diagnosed_oldfemales.ci_low = diagnosed_old_females.ci_low / popsize_old_females * 1e5,
    diagnosed_youngmales.ci_low = diagnosed_young_males.ci_low / popsize_young_males * 1e5,
    diagnosed_oldmales.ci_low = diagnosed_old_males.ci_low / popsize_old_males * 1e5,

    diagnosed_youngfemales.median = diagnosed_young_females.median / popsize_young_females * 1e5,
    diagnosed_oldfemales.median = diagnosed_old_females.median / popsize_old_females * 1e5,
    diagnosed_youngmales.median = diagnosed_young_males.median / popsize_young_males * 1e5,
    diagnosed_oldmales.median = diagnosed_old_males.median / popsize_old_males * 1e5) %>% 

  # remove population sizes
  select(
    -c(popsize_young_females, 
       popsize_young_males,
       popsize_old_females, 
       popsize_old_males,
       diagnosed_young_females.median,
       diagnosed_old_females.median,
       diagnosed_young_males.median,
       diagnosed_old_males.median,
       incidence_young_females.median,
       incidence_old_females.median,
       incidence_young_males.median,
       incidence_old_males.median,
       prevalence_young_females.median,
       prevalence_old_females.median,
       prevalence_young_males.median,
       prevalence_old_males.median,
       diagnosed_young_females.ci_high,
       diagnosed_old_females.ci_high,
       diagnosed_young_males.ci_high,
       diagnosed_old_males.ci_high,
       incidence_young_females.ci_high,
       incidence_old_females.ci_high,
       incidence_young_males.ci_high,
       incidence_old_males.ci_high,
       prevalence_young_females.ci_high,
       prevalence_old_females.ci_high,
       prevalence_young_males.ci_high,
       prevalence_old_males.ci_high,
       diagnosed_young_females.ci_low,
       diagnosed_old_females.ci_low,
       diagnosed_young_males.ci_low,
       diagnosed_old_males.ci_low,
       incidence_young_females.ci_low,
       incidence_old_females.ci_low,
       incidence_young_males.ci_low,
       incidence_old_males.ci_low,
       prevalence_young_females.ci_low,
       prevalence_old_females.ci_low,
       prevalence_young_males.ci_low,
       prevalence_old_males.ci_low
       )) %>% 

  # reshape the columns into rows with column-variables
  reshape2::melt(id.vars = c('state', 'intervention', 'year'), value.name = 'rate') -> df2_reshaped

  # separate the variable column into variable and sex
  df2_reshaped %>% tidyr::separate(col = 'variable', into = c('group', 'statistic'), sep = '\\.', fill = 'left') -> df2_reshaped_1

  df2_reshaped_1 %>% tidyr::separate(col = 'group', into = c('variable', 'group'), extra='drop',
    sep='_') -> df3

  df3 %>% tidyr::spread(statistic, rate) -> df3_spread

  saveRDS(df3_spread, file.path(system.file('interventions/', package='syphilis'), 'intvs_formatted.rds'))

  
  # ggplot(filter(df3_spread, state == 'Louisiana', variable == 'prevalence', group == 'oldmales'),
  #   aes(x = year, y = mean, ymax = ci_high, ymin = ci_low, fill = intervention, color = intervention)) + 
  #   geom_ribbon(size = 0, alpha = 0.8) + geom_line() + 
  #   theme_minimal() + 
  #   facet_wrap(~intervention)


}


compute_interventions_change_in_outcomes <- function() { 

  # read in raw data
  la_df <- readRDS(system.file('interventions/la_interventions_df.rds', package='syphilis'))
  ma_df <- readRDS(system.file('interventions/ma_interventions_df.rds', package='syphilis'))

  # > colnames(la_df)
  #  [1] "sim"                      "intervention"             "year"                     "popsize_all"              "popsize_young_females"    "popsize_old_females"      "popsize_young_males"
  #  [8] "popsize_old_males"        "diagnosed_all"            "diagnosed_young_females"  "diagnosed_old_females"    "diagnosed_young_males"    "diagnosed_old_males"      "incidence_all"
  # [15] "incidence_young_females"  "incidence_old_females"    "incidence_young_males"    "incidence_old_males"      "prevalence_all"           "prevalence_young_females" "prevalence_old_females"
  # [22] "prevalence_young_males"   "prevalence_old_males"

  # combine the data frames
  la_df <- cbind.data.frame(state = 'Louisiana', la_df)
  ma_df <- cbind.data.frame(state = 'Massachusetts', ma_df)
  df <- rbind.data.frame(la_df, ma_df)

  # we're only going to compute change from basecase in the last intervention year
  df %<>% filter(year == max(year))

  # get just basecase
  basecase_df <- filter(df, intervention == 'basecase')

  # filter out basecase
  df <- filter(df, intervention != 'basecase')

  # put intervention data and basecase data next to each other for the same sim, year, state
  df2 <- merge(df, basecase_df, by = c('state', 'sim', 'year'), suffixes = c('', '.basecase'))

  # intervention.basecase doesn't contain any information, just 'basecase'
  df2 %<>% select(-intervention.basecase)

  # popsize is constant, so it cancels out of the change compared to basecase equation
  df2 %<>% select(-starts_with('popsize'))

  # list of variables to compare
  varlist <- c("incidence_all", "prevalence_all", "diagnosed_all")

  for (var in varlist) {  # compute change from basecase
    df2[[var]] <- (df2[[var]] - df2[[paste0(var, '.basecase')]]) / df2[[paste0(var, '.basecase')]]
  }

  df2 %<>% select(-ends_with('.basecase')) # drop basecase values

  df2[,varlist] <- df2[,varlist]*100 # format as a percent change

  # save raw data
  saveRDS(df2, file.path(system.file("interventions/", package='syphilis'), 'intvs_change_from_basecase.rds'))

  # function to compute sixnums (mean, median, interquartile, confidence intervals) for a list of columns
  summarize_sixnum <- function(df, cols) { 
    df %>% 
      group_by(state, year, intervention) %>% 
      summarize_at(
        vars(cols), list(
          mean = ~mean(., na.rm=T),
          median = ~quantile(., .5, na.rm=T),
          ci_high = ~quantile(., .975, na.rm=T),
          ci_low = ~quantile(., .025, na.rm=T),
          iq_high = ~quantile(., .75, na.rm=T),
          iq_low = ~quantile(., .25, na.rm=T)
          ))
  }

  # compute sixnums
  df3 <- summarize_sixnum(df2, varlist)

  # replace _ with . for easier splitting later
  colnames(df3) %<>%
    gsub(pattern = "_mean", replacement = ".mean", .) %>% 
    gsub(pattern = "_median", replacement = ".median", .) %>% 
    gsub(pattern = "_ci_high", replacement = ".ci_high", .) %>% 
    gsub(pattern = "_ci_low", replacement = ".ci_low", .) %>% 
    gsub(pattern = "_iq_high", replacement = ".iq_high", .) %>% 
    gsub(pattern = "_iq_low", replacement = ".iq_low", .)

  df4 <- df3 %>% tidyr::gather(key = "outcome", value = "change", -state, -year, -intervention)

  df4$outcome %<>% 
    gsub(pattern = "_all", replacement = "", .)

  df4 %<>% tidyr::separate(col = 'outcome', into = c('outcome', 'statistic'), sep = '\\.')

  df4 %<>% tidyr::spread('statistic', 'change')

  saveRDS(df4, file.path(system.file("interventions/", package='syphilis'), 
    'intvs_change_from_basecase_summarized.rds'))
    
}



simulate_outcomes_by_sex_by_state <- function(state) { 
  
  state <<- state
  # load state model
  load.start()

  # get state parametrization
  theta <- get(paste0('theta_', tolower(state)))

  # simulate interventions, record diagnoses by sex
  df <- simulate_interventions(theta, info_extractor = extract_outcomes)

  df %<>% 
  # reshape data into rates per 100k
  dplyr::mutate(
    # prevalence 
    # prevalence_all = prevalence_all / popsize_all * 1e5,
    # prevalence_young = prevalence_young / popsize_young * 1e5,
    # prevalence_old = prevalence_old / popsize_old * 1e5,

    # prevalence_females = prevalence_females / popsize_females * 1e5,
    prevalence_youngfemales = prevalence_young_females / popsize_young_females * 1e5,
    prevalence_oldfemales = prevalence_old_females / popsize_old_females * 1e5,

    # prevalence_males = prevalence_males / popsize_males * 1e5,
    prevalence_youngmales = prevalence_young_males / popsize_young_males * 1e5,
    prevalence_oldmales = prevalence_old_males / popsize_old_males * 1e5,

    # prevalence_msw = prevalence_msw / popsize_msw * 1e5,
    # prevalence_youngmsw = prevalence_young_msw / popsize_young_msw * 1e5,
    # prevalence_oldmsw = prevalence_old_msw / popsize_old_msw * 1e5,

    # prevalence_msm = prevalence_msm / popsize_msm * 1e5,
    # prevalence_youngmsm = prevalence_young_msm / popsize_young_msm * 1e5,
    # prevalence_oldmsm = prevalence_old_msm / popsize_old_msm * 1e5,

    # prevalence_new = prevalence_new / popsize_all * 1e5,
    # prevalence_reinfected = prevalence_reinfected / popsize_all * 1e5,
    # prevalence_high_activity = prevalence_high_activity / popsize_all * 1e5,
    # prevalence_low_activity = prevalence_low_activity / popsize_all * 1e5,

    # incidence 
    # incidence_all = incidence_all / popsize_all * 1e5,
    # incidence_young = incidence_young / popsize_young * 1e5,
    # incidence_old = incidence_old / popsize_old * 1e5,

    # incidence_females = incidence_females / popsize_females * 1e5,
    incidence_youngfemales = incidence_young_females / popsize_young_females * 1e5,
    incidence_oldfemales = incidence_old_females / popsize_old_females * 1e5,

    # incidence_males = incidence_males / popsize_males * 1e5,
    incidence_youngmales = incidence_young_males / popsize_young_males * 1e5,
    incidence_oldmales = incidence_old_males / popsize_old_males * 1e5,

    # incidence_msw = incidence_msw / popsize_msw * 1e5,
    # incidence_youngmsw = incidence_young_msw / popsize_young_msw * 1e5,
    # incidence_oldmsw = incidence_old_msw / popsize_old_msw * 1e5,

    # incidence_msm = incidence_msm / popsize_msm * 1e5,
    # incidence_youngmsm = incidence_young_msm / popsize_young_msm * 1e5,
    # incidence_oldmsm = incidence_old_msm / popsize_old_msm * 1e5,

    # incidence_reinfected = incidence_reinfected / popsize_all * 1e5,
    # incidence_new = incidence_new / popsize_all * 1e5,
    # incidence_high_activity = incidence_high_activity / popsize_all * 1e5,
    # incidence_low_activity = incidence_low_activity / popsize_all * 1e5,

    # diagnoses
    # diagnosed_all = diagnosed_all / popsize_all * 1e5,
    # diagnosed_young = diagnosed_young / popsize_young * 1e5,
    # diagnosed_old = diagnosed_old / popsize_old * 1e5,

    # diagnosed_females = diagnosed_females / popsize_females * 1e5,
    diagnosed_youngfemales = diagnosed_young_females / popsize_young_females * 1e5,
    diagnosed_oldfemales = diagnosed_old_females / popsize_old_females * 1e5,

    # diagnosed_males = diagnosed_males / popsize_males * 1e5,
    diagnosed_youngmales = diagnosed_young_males / popsize_young_males * 1e5,
    diagnosed_oldmales = diagnosed_old_males / popsize_old_males * 1e5) %>% 

    # diagnosed_msw = diagnosed_msw / popsize_msw * 1e5,
    # diagnosed_youngmsw = diagnosed_young_msw / popsize_young_msw * 1e5,
    # diagnosed_oldmsw = diagnosed_old_msw / popsize_old_msw * 1e5,

    # diagnosed_msm = diagnosed_msm / popsize_msm * 1e5,
    # diagnosed_youngmsm = diagnosed_young_msm / popsize_young_msm * 1e5,
    # diagnosed_oldmsm = diagnosed_old_msm / popsize_old_msm * 1e5) %>% 

  # remove population sizes
  select(
    -c(popsize_young_females, 
       popsize_young_males,
       popsize_old_females, 
       popsize_old_males
       )) %>% 

  # reshape the columns into rows with column-variables
  reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'rate') %>% 

  # separate the variable column into variable and sex
  tidyr::separate(col = variable, into = c('variable', 'group'), extra='drop')

  return(df)
}


# simulate_diagnoses_by_sex_by_state <- function(state) { 
  
#   state <<- state
#   # load state model
#   load.start()

#   # get state parametrization
#   theta <- get(paste0('theta_', tolower(state)))

#   # simulate interventions, record diagnoses by sex
#   df <- simulate_interventions(theta, info_extractor = extract_diagnoses_by_sex)

#   # reshape data into diagnoses per 100k
#   df %<>% 
#   mutate(diagnosed_females_per_100k = diagnosed_females / popsize_females * 100000,
#          diagnosed_msw_per_100k = diagnosed_msw / popsize_msw * 100000,
#          diagnosed_msm_per_100k = diagnosed_msm / popsize_msm * 100000) %>% 
#   select(-c(diagnosed_females, diagnosed_msw, diagnosed_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
#   reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'diagnosis_rate') %>% 
#   tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
#   select(-prefix)

#   return(df)
# }


# simulate_incidence_by_sex_by_state <- function(state) { 
  
#   state <<- state
#   # load state model
#   load.start()

#   # get state parametrization
#   theta <- get(paste0('theta_', tolower(state)))

#   # simulate interventions, record diagnoses by sex
#   df <- simulate_interventions(theta, info_extractor = extract_incidence_by_sex)

#   # reshape data into diagnoses per 100k
#   df %<>% 
#   mutate(incidence_females_per_100k = incidence_females / popsize_females * 100000,
#          incidence_msw_per_100k = incidence_msw / popsize_msw * 100000,
#          incidence_msm_per_100k = incidence_msm / popsize_msm * 100000) %>% 
#   select(-c(incidence_females, incidence_msw, incidence_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
#   reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'incidence_rate') %>% 
#   tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
#   select(-prefix)

#   return(df)
# }


# simulate_prevalence_by_sex_by_state <- function(state) { 
  
#   state <<- state
#   # load state model
#   load.start(model.end = 115)

#   # get state parametrization
#   theta <- get(paste0('theta_', tolower(state)))

#   # simulate interventions, record diagnoses by sex
#   df <- simulate_interventions(theta, info_extractor = extract_prevalence_by_sex)

#   # reshape data into diagnoses per 100k
#   df %<>% 
#   mutate(prevalence_females_per_100k = prevalence_females / popsize_females * 100000,
#          prevalence_msw_per_100k = prevalence_msw / popsize_msw * 100000,
#          prevalence_msm_per_100k = prevalence_msm / popsize_msm * 100000) %>% 
#   select(-c(prevalence_females, prevalence_msw, prevalence_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
#   reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'prevalence_rate') %>% 
#   tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
#   select(-prefix)

#   return(df)
# }

simulate_outcomes_by_sex_both_states <- function() { 
  rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', simulate_outcomes_by_sex_by_state('LA')),
    cbind.data.frame(state = 'Massachusetts', simulate_outcomes_by_sex_by_state('MA')))
}

plot_state_outcomes_by_sex_and_age <- function(state, outcome, years, include_legend = FALSE, include_ci = FALSE, scales = 'fixed') { 

  if (missing(years)) {
    min_year <- intervention_start_gregorian
    max_year <- 2017 # intervention_stop_gregorian
  }

  if (! exists('df') || 
  colnames(df) != c("state", "intervention", "year", "variable", "group", "ci_high",
"ci_low", "median")) { 
    stop("df should be read in from inst/interventions/intvs_formatted.rds")
  }

  df %<>% filter(year <= max_year, year >= min_year)

  state1 <- enquo(state)
  outcome1 <- enquo(outcome)

  # group_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM', males = 'M', all = 'All', 
  #   new = 'New Infections', reinfected = 'Reinfected', young = 'Young', youngfemales = 'Young F', 
  #   youngmales = 'Young M', youngmsw = 'Young MSW', youngmsm = 'Young MSM',
  #   old = 'Old', oldmales = 'Old M', oldfemales = 'Old F', oldmsm = 'Old MSM', oldmsw = 'Old MSW')

  outcome_lookup <- c(prevalence = "Prevalent early syphilis infections\nper 100,000 population\n", 
    incidence = "Incident syphilis cases\nper 100,000 population\n",
    diagnosed = "Reported early syphilis cases\nper 100,000 population\n")

  df %<>% mutate(group = recode(group, youngmales = "M: 25-44 y", youngfemales = "F: 25-44 y", oldmales = "M: 45-64 y", oldfemales = "F: 45-64 y"))

  groupnames <- c("M: 25-44 y", "F: 25-44 y", "M: 45-64 y", "F: 45-64 y")

  ggplot(
    data = filter(df, state == !! state1, group %in% groupnames, variable == !! outcome1), 
    mapping = aes(x = year, y = median, ymax = ci_high, ymin = ci_low, color = intervention, fill = intervention, linetype = intervention)) + 
    geom_line(size = 1, alpha=0.8) + 
    # if (include_ci) { geom_ribbon(size = 0, alpha=0.5) } else { NULL } + 
    ylab('') + 
    xlab('') + 
    expand_limits(y=0) + 
    facet_wrap(group~., scales=scales) + 
    theme_minimal() + 
    scale_color_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        # 'Prior Diagnosis Annual', 
        'Prior Diagnosis every 3 months', 
        # 'High Activity Annual',
        'High Activity every 3 months'), 
        # values = c("#05141E", "#762B19", "#3D507A", "#657062", "#D14E3E", "#E78A40", "#EBD799")
        # values = c("#7BA46C", "#602D31", "#008D91", "#0A789F", "#C6A28A", "#61B8D3", "#EACF9E")
        # values = c("#7FC97F", "#BEAED4", "#FDC086", "#AAAA99", "#386CB0", "#F0027F", "#BF5B17")
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_fill_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        # 'Prior Diagnosis Annual', 
        'Prior Diagnosis every 3 months', 
        # 'High Activity Annual',
        'High Activity every 3 months'), 
        # values = c("#05141E", "#762B19", "#3D507A", "#657062", "#D14E3E", "#E78A40", "#EBD799")
        # values = c("#7BA46C", "#602D31", "#008D91", "#0A789F", "#C6A28A", "#61B8D3", "#EACF9E")
        # values = c("#7FC97F", "#BEAED4", "#FDC086", "#AAAA99", "#386CB0", "#F0027F", "#BF5B17")
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
                   

      scale_linetype_manual(
        labels = c('Base Case', 
          'Guidelines in MSM', 
          'MSM every 3 months', 
          'Prior Diagnosis every 3 months', 
          'High Activity every 3 months'), 
        values = c(1,5,2,3,4))  +
      # scale_x_continuous(breaks = seq(intervention_start_gregorian,
      #   intervention_stop_gregorian, by=2)) + 
			labs( x="Year", y=outcome_lookup[[outcome]], title=state, color = 'Intervention', linetype = 'Intervention', fill = 'Intervention') +
    theme(legend.position= if(include_legend) "right" else "none",  axis.text.x=element_text(angle=45,hjust=0.5,size=8), axis.text.y=element_text(size=10),  
      axis.title.x = element_text(size = 10),
			plot.title=element_text(size=12, hjust=0.5),
			text=element_text(size=12)) + 
    if (include_ci) { geom_ribbon(size = 0, alpha = 0.5) } 
      # if (include_title) { 
      #   ggtitle(paste0(outcome_lookup[[outcome]], ' per 100,000'), subtitle = paste0(group_lookup[[group]], " - ", state)) 
      # } else { NULL } 
}


plot_state_outcomes_by_sex_and_age_conf_ints <- function(state, outcome, years, group, include_legend = FALSE) { 

  if (missing(years)) {
    min_year <- intervention_start_gregorian
    max_year <- 2017 # intervention_stop_gregorian
  }

  if (! exists('df') || 
  colnames(df) != c(
  "state","intervention","year","variable","group","ci_high","ci_low","mean")) { 
    stop("df should be read in from inst/interventions/intvs_formatted.rds")
  }

  df %<>% filter(year <= max_year, year >= min_year)

  state1 <- enquo(state)
  outcome1 <- enquo(outcome)
  group1 <- enquo(group)

  # group_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM', males = 'M', all = 'All', 
  #   new = 'New Infections', reinfected = 'Reinfected', young = 'Young', youngfemales = 'Young F', 
  #   youngmales = 'Young M', youngmsw = 'Young MSW', youngmsm = 'Young MSM',
  #   old = 'Old', oldmales = 'Old M', oldfemales = 'Old F', oldmsm = 'Old MSM', oldmsw = 'Old MSW')

  outcome_lookup <- c(prevalence = "Prevalent early syphilis infections\nper 100,000 population\n", 
    incidence = "Incident syphilis cases\nper 100,000 population\n",
    diagnosed = "Reported early syphilis cases\nper 100,000 population\n")

  df %<>% mutate(group = recode(group, youngmales = "M: 25-44 y", youngfemales = "F: 25-44 y", oldmales = "M: 45-64 y", oldfemales = "F: 45-64 y"))

  df %<>% mutate(intervention = recode(intervention, 
    basecase = 'Base Case', msm_annual_hr_msm_quarterly = 'Guidelines in MSM',
    msm_quarterly = 'MSM every 3 months', prior_diagnosis_quarterly = 'Prior Diagnosis every 3 months',
    high_activity_quarterly = 'High Activity every 3 months'))

  groupnames <- c("M: 25-44 y", "F: 25-44 y", "M: 45-64 y", "F: 45-64 y")
  group_lookup <- c(youngfemales = 'F: 25-44 y', youngmales = 'M: 25-44 y', oldfemales = 'F: 45-64 y', oldmales = 'M: 45-64 y')

  ggplot(
    data = filter(df, state == !! state1, group == group_lookup[[!! group1]], variable == !! outcome1), 
    mapping = aes(x = year, y = mean, ymax = ci_high, ymin = ci_low, fill =
      intervention, color = intervention, linetype = intervention)) + 

    geom_line(size = 1, alpha=0.8) + 
    geom_ribbon(size = 0, alpha=0.5) + 

    ylab('') + 
    xlab('') + 
    expand_limits(y=c(0, max(df$ci_high, na.rm=T))) + 
    facet_wrap(intervention~.) + 
    theme_minimal() + 
    scale_color_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
                   
    scale_fill_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 

    scale_linetype_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months', 
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
      values = c(1,5,2,3,4))  +

    labs( x="Year", y=outcome_lookup[[outcome]], title=state, subtitle = group_lookup[[group]], 
      color = 'Intervention', linetype = 'Intervention') +

    theme(legend.position= if(include_legend) "right" else "none",
      axis.text.x=element_text(angle=45,hjust=0.5,size=8),
      axis.text.y=element_text(size=10),  
      axis.title.x = element_text(size = 10),
      plot.title=element_text(size=12, hjust=0.5),
      text=element_text(size=12))
}




#' @import tidyverse
# simulate_and_plot_diagnoses_by_sex_both_states <- function() { 
# 	load_globals()

#   la_df <- simulate_diagnoses_by_sex_by_state('LA')
#   ma_df <- simulate_diagnoses_by_sex_by_state('MA')

#   df <- rbind.data.frame(
#     cbind.data.frame(state = 'Louisiana', la_df),
#     cbind.data.frame(state = 'Massachusetts', ma_df))

#   df %<>% filter(year <= 2020, year >= 2012)

#   produce_plot <- function(state, sex) { 
#     state1 <- enquo(state)
#     sex1 <- enquo(sex)

#     sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

#     ggplot( data = filter(df, state == !! state, group == !! sex), 
#             mapping = aes(x = year, y = diagnosis_rate, color = intervention, linetype = intervention)) + 
#       geom_line(size = 1) + 
#       ylab('') + 
#       xlab('') + 
#       expand_limits(y=0) + 
#       scale_color_manual(values = c('red', 'forestgreen', 'skyblue', 'orange', 'coral')) + 
#       scale_x_continuous(breaks = seq(2012, 2020, by=2)) + 
#       theme_bw() + 
#       ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
#   }


#   p1 <- produce_plot('Louisiana', 'females') 
#   p2 <- produce_plot('Massachusetts', 'females')
#   p3 <- produce_plot('Louisiana', 'msw')
#   p4 <- produce_plot('Massachusetts', 'msw')
#   p5 <- produce_plot('Louisiana', 'msm')
#   p6 <- produce_plot('Massachusetts', 'msm')

#   legend <- cowplot::get_legend(p1)

#   p1 <- p1 + theme(legend.position = 'none')
#   p2 <- p2 + theme(legend.position = 'none')
#   p3 <- p3 + theme(legend.position = 'none')
#   p4 <- p4 + theme(legend.position = 'none')
#   p5 <- p5 + theme(legend.position = 'none')
#   p6 <- p6 + theme(legend.position = 'none')

#   plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

#   y.grob <- textGrob("Diagnosed Cases Per 100,000", 
#                    gp=gpar(fontface="bold", fontsize=15), rot=90)

#   x.grob <- textGrob("Year", 
#                    gp=gpar(fontface="bold", fontsize=15))

#   plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

#   # ggsave(plot = plt, 
#   #   filename = file.path(system.file("figures/diagnoses_by_sex_in_interventions/", package = 'syphilis'), 'diagnoses_by_sex_in_interventions.png'),
#   #   width = 10)

#   return(plt)
# }


#' @import tidyverse
# simulate_and_plot_incidence_by_sex_both_states <- function() { 
# 	load_globals(model.end = 116)

#   la_df <- simulate_incidence_by_sex_by_state('LA')
#   ma_df <- simulate_incidence_by_sex_by_state('MA')

#   df <- rbind.data.frame(
#     cbind.data.frame(state = 'Louisiana', la_df),
#     cbind.data.frame(state = 'Massachusetts', ma_df))

#   df %<>% filter(year <= 2020, year >= 2012)

#   produce_plot <- function(state, sex) { 
#     state1 <- enquo(state)
#     sex1 <- enquo(sex)

#     sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

#     ggplot( data = filter(df, state == !! state, group == !! sex), 
#             mapping = aes(x = year, y = incidence_rate, color = intervention, linetype = intervention)) + 
#       geom_line(size = 1) + 
#       ylab('') + 
#       xlab('') + 
#       expand_limits(y=0) + 
#       scale_color_manual(
#         labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'), 
#         values = c('red', 'forestgreen', 'skyblue', 'orange', 'purple')) + 
#       scale_linetype_manual(
#         labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'),
#         values = c(1,5,2,3,4))  +
#       scale_x_continuous(breaks = seq(2012, 2020, by=2)) + 
#       theme_bw() + 
      
#       ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
#   }


#   p1 <- produce_plot('Louisiana', 'females')
#   p2 <- produce_plot('Massachusetts', 'females')
#   p3 <- produce_plot('Louisiana', 'msw')
#   p4 <- produce_plot('Massachusetts', 'msw')
#   p5 <- produce_plot('Louisiana', 'msm')
#   p6 <- produce_plot('Massachusetts', 'msm')

#   legend <- cowplot::get_legend(p1)

#   p1 <- p1 + theme(legend.position = 'none')
#   p2 <- p2 + theme(legend.position = 'none')
#   p3 <- p3 + theme(legend.position = 'none')
#   p4 <- p4 + theme(legend.position = 'none')
#   p5 <- p5 + theme(legend.position = 'none')
#   p6 <- p6 + theme(legend.position = 'none')

#   plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

#   y.grob <- textGrob("Incidence Per 100,000", 
#                    gp=gpar(fontface="bold", fontsize=15), rot=90)

#   x.grob <- textGrob("Year", 
#                    gp=gpar(fontface="bold", fontsize=15))

#   plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

#   ggsave(plot = plt, 
#     filename = file.path(system.file("figures/", package = 'syphilis'), 'incidence_by_sex_in_interventions.png'),
#     width = 10)

#   return(plt)
# }

#' @import tidyverse
# simulate_and_plot_prevalence_by_sex_both_states <- function() { 
# 	load_globals(model.end = 110)

#   la_df <- simulate_prevalence_by_sex_by_state('LA')
#   ma_df <- simulate_prevalence_by_sex_by_state('MA')

#   df <- rbind.data.frame(
#     cbind.data.frame(state = 'Louisiana', la_df),
#     cbind.data.frame(state = 'Massachusetts', ma_df))

#   df %<>% filter(year <= 2020, year >= 2012)

#   produce_plot <- function(state, sex) { 
#     state1 <- enquo(state)
#     sex1 <- enquo(sex)

#     sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

#     ggplot( data = filter(df, state == !! state, group == !! sex), 
#             mapping = aes(x = year, y = prevalence_rate, color = intervention, linetype = intervention)) + 
#       geom_line(size = 1) + 
#       ylab('') + 
#       xlab('') + 
#       expand_limits(y=0) + 
#       scale_color_manual(
#         labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'), 
#         values = c('red', 'forestgreen', 'skyblue', 'orange', 'purple')) + 
#       scale_linetype_manual(
#         labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'),
#         values = c(1,5,2,3,4))  +
#       scale_x_continuous(breaks = seq(2012, 2020, by=2)) + 
#       theme_bw() + 
      
#       ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
#   }


#   p1 <- produce_plot('Louisiana', 'females')
#   p2 <- produce_plot('Massachusetts', 'females')
#   p3 <- produce_plot('Louisiana', 'msw')
#   p4 <- produce_plot('Massachusetts', 'msw')
#   p5 <- produce_plot('Louisiana', 'msm')
#   p6 <- produce_plot('Massachusetts', 'msm')

#   legend <- cowplot::get_legend(p1)

#   p1 <- p1 + theme(legend.position = 'none')
#   p2 <- p2 + theme(legend.position = 'none')
#   p3 <- p3 + theme(legend.position = 'none')
#   p4 <- p4 + theme(legend.position = 'none')
#   p5 <- p5 + theme(legend.position = 'none')
#   p6 <- p6 + theme(legend.position = 'none')

#   plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

#   y.grob <- textGrob("Prevalence Per 100,000", 
#                    gp=gpar(fontface="bold", fontsize=15), rot=90)

#   x.grob <- textGrob("Year", 
#                    gp=gpar(fontface="bold", fontsize=15))

#   plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

#   ggsave(plot = plt, 
#     filename = file.path(system.file("figures/", package = 'syphilis'), 'prevalence_by_sex_in_interventions.png'),
#     width = 10)

#   return(plt)
# }





breakdown_basecase <- function() { 
  
  for (state in c('LA', 'MA')) { 
    theta <- get(paste0('theta_', tolower(state)))
    state <<- state
    load.start()

		e <- constructSimulationEnvironment(theta)
		e$params$output_weekly <- TRUE
		runSimulation(e)
		assign(paste0(tolower(state), '_df'), extract_outcomes(e$sim$out))
  }

  df <- rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', la_df),
    cbind.data.frame(state = 'Massachusetts', ma_df))

  df %<>% 
    mutate(
    # prevalence 
    prevalence_all = prevalence_all / popsize_all * 1e5,
    prevalence_females = prevalence_females / popsize_all * 1e5,
    prevalence_msw = prevalence_msw / popsize_all * 1e5,
    prevalence_msm = prevalence_msm / popsize_all * 1e5,
    prevalence_new = prevalence_new / popsize_all * 1e5,
    prevalence_reinfected = prevalence_reinfected / popsize_all * 1e5,
    prevalence_high_activity = prevalence_high_activity / popsize_all * 1e5,
    prevalence_low_activity = prevalence_low_activity / popsize_all * 1e5,
    prevalence_noninfectious = prevalence_noninfectious / popsize_all * 1e5,
    prevalence_infectious = prevalence_infectious / popsize_all * 1e5,

    # incidence 
    incidence_all = incidence_all / popsize_all * 1e5,
    incidence_females = incidence_females / popsize_all * 1e5,
    incidence_msw = incidence_msw / popsize_all * 1e5,
    incidence_msm = incidence_msm / popsize_all * 1e5,
    incidence_reinfected = incidence_reinfected / popsize_all * 1e5,
    incidence_new = incidence_new / popsize_all * 1e5,
    incidence_high_activity = incidence_high_activity / popsize_all * 1e5,
    incidence_low_activity = incidence_low_activity / popsize_all * 1e5,

    # diagnoses
    diagnosed_all = diagnosed_all / popsize_all * 1e5,
    diagnosed_females = diagnosed_females / popsize_all * 1e5,
    diagnosed_msw = diagnosed_msw / popsize_all * 1e5,
    diagnosed_msm = diagnosed_msm / popsize_all * 1e5) %>% 

  # remove population sizes
  select(
    -c(popsize_females, 
       popsize_msw, 
       popsize_msm,
       tests_all)) %>% 

  # reshape the columns into rows with column-variables
  reshape2::melt(id.vars = c('state', 'year'), value.name = 'rate') %>% 

  # separate the variable column into variable and sex
  tidyr::separate(col = variable, into = c('variable', 'group'), extra='drop')

  return(df)

}


# make_intervention_analysis_plots <- function() { 
#   df <- simulate_outcomes_by_sex_both_states()

#   # df %<>% filter((intervention == 'basecase' & year %% 1 == 0 | year < 2012) | intervention != 'basecase')

#   la_prevalence <- plot_state_outcomes_by_sex('Louisiana', 'all', 'prevalence')

#   ma_prevalence <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'prevalence')

#   ggsave(plot=la_prevalence, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence.png'),
#     width = 10, height = 6)


#   ggsave(plot=ma_prevalence, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_prevalence.png'),
#     width = 10, height = 6)


#   la_incidence <- plot_state_outcomes_by_sex('Louisiana', 'all', 'incidence')
#   ma_incidence <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'incidence')

#   ggsave(plot=la_incidence, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence.png'),
#     width = 10, height = 6)

#   ggsave(plot=ma_incidence, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence.png'),
#     width = 10, height = 6)

#   la_diagnosed <- plot_state_outcomes_by_sex('Louisiana', 'all', 'diagnosed')
#   ma_diagnosed <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'diagnosed')

#   ggsave(plot=la_diagnosed, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_diagnosed.png'),
#     width = 10, height = 6)

#   ggsave(plot=ma_diagnosed, 
#     file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_diagnosed.png'),
#     width = 10, height = 6)

#   df <- breakdown_basecase()

#   plot_breakdown_by_state_outcome_grouping <- function(state, outcome, grouping) { 
#     state1 <- enquo(state)
#     outcome1 <- enquo(outcome)
#     grouping1 <- enquo(grouping)

#     grouping_lookup <- list('risk' = c('high', 'low'), sex = c('msm', 'msw', 'females'), reinfection = c('new', 'reinfected'), infectiousness = c('noninfectious', 'infectious'))
    

#     grouping_color_lookup <- c(risk = 'Reds', reinfection = 'Purples', sex = 'Blues', infectiousness = 'Oranges')

#     ggplot(
#         data = filter(df, year >= 2012, year <= 2016, state == !! state1, variable == !! outcome1, group %in% grouping_lookup[[grouping]]),
#         mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#         geom_area(alpha=0.8) + 
#         scale_fill_brewer(palette = grouping_color_lookup[[grouping]]) + 
#         theme_bw() + 
#         scale_x_continuous(breaks = seq(2012, 2016, by=2)) + 
#         ggtitle(paste0(tools::toTitleCase(outcome), " per 100,000 in ", state), subtitle = paste0("Broken down by ", tools::toTitleCase(grouping)))
#   }

#   la_prevalence_breakdown_risk <- ggplot(
#     data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('high', 'low')),
#     mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#     geom_area(alpha=0.5) + 
#     scale_fill_brewer(palette = 'Reds') + 
#     theme_bw()

#   la_prevalence_breakdown_infectiousness <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'prevalence', 'infectiousness')
#   ggsave(plot = la_prevalence_breakdown_infectiousness, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence_breakdown_infectiousness.png'))

#   ma_prevalence_breakdown_infectiousness <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'prevalence', 'infectiousness')
#   ggsave(plot = ma_prevalence_breakdown_infectiousness, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_prevalence_breakdown_infectiousness.png'))

#   la_incidence_breakdown_risk <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'risk')
#   ggsave(plot = la_incidence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_risk.png'))

#   ma_incidence_breakdown_risk <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'risk')
#   ggsave(plot = ma_incidence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_risk.png'))

#   la_incidence_breakdown_sex <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'sex')
#   ggsave(plot = la_incidence_breakdown_sex, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_sex.png'))

#   ma_incidence_breakdown_sex <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'sex')
#   ggsave(plot = ma_incidence_breakdown_sex, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_sex.png'))


#   la_incidence_breakdown_reinfection <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'reinfection')
#   ggsave(plot = la_incidence_breakdown_reinfection, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_reinfection.png'))

#   ma_incidence_breakdown_reinfection <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'reinfection')
#   ggsave(plot = ma_incidence_breakdown_reinfection, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_reinfection.png'))
#   # la_prevalence_breakdown_sex <- ggplot(
#   #   data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('msm', 'msw', 'females')),
#   #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#   #   geom_area(alpha=0.5) + 
#   #   scale_fill_brewer(palette = 'Reds') + 
#   #   theme_bw()

#   # la_prevalence_breakdown_reinfection <- ggplot(
#   #   data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('new', 'reinfected')),
#   #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#   #   geom_area(alpha=0.5) + 
#   #   scale_fill_brewer(palette = 'Purples') + 
#   #   theme_bw()

#   # ma_prevalence_breakdown_risk <- ggplot(
#   #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('high', 'low')),
#   #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#   #   geom_area(alpha=0.5) + 
#   #   scale_fill_brewer(palette = 'Reds') + 
#   #   theme_bw()

#   # ma_prevalence_breakdown_sex <- ggplot(
#   #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('msm', 'msw', 'females')),
#   #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#   #   geom_area(alpha=0.5) + 
#   #   scale_fill_brewer(palette = 'Reds') + 
#   #   theme_bw()

#   # ma_prevalence_breakdown_reinfection <- ggplot(
#   #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('new', 'reinfected')),
#   #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
#   #   geom_area(alpha=0.5) + 
#   #   scale_fill_brewer(palette = 'Purples') + 
#   #   theme_bw()


#   ggsave(plot = la_prevalence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence_breakdown_risk.png'))


# }


# make_4x2_young_group_by_state_outcomes_plot <- function(outcome) { 

#   require(cowplot)

#   plots <- list()

#   plots[['la_yf']] <- plot_state_outcomes_by_sex('Louisiana', 'youngfemales', outcome)
#   plots[['la_ym']] <- plot_state_outcomes_by_sex('Louisiana', 'youngmales', outcome)
#   plots[['la_ymsm']] <- plot_state_outcomes_by_sex('Louisiana', 'youngmsm', outcome)
#   plots[['la_ymsw']] <- plot_state_outcomes_by_sex('Louisiana', 'youngmsw', outcome)

#   plots[['ma_yf']] <- plot_state_outcomes_by_sex('Massachusetts', 'youngfemales', outcome)
#   plots[['ma_ym']] <- plot_state_outcomes_by_sex('Massachusetts', 'youngmales', outcome)
#   plots[['ma_ymsm']] <- plot_state_outcomes_by_sex('Massachusetts', 'youngmsm', outcome)
#   plots[['ma_ymsw']] <- plot_state_outcomes_by_sex('Massachusetts', 'youngmsw', outcome)

#   titles <- c('Young F', 
# 							'Young M', 
# 							'Young MSM', 
# 							'Young MSW', 
# 							'Young F', 
# 							'Young M', 
# 							'Young MSM', 
# 							'Young MSW')

# 	la_title <- ggdraw() + 
#     draw_label(
# 			"Louisiana",
# 			fontface = 'bold', hjust = 0.5)

# 	ma_title <- ggdraw() + 
#     draw_label(
# 			"Massachusetts",
# 			fontface = 'bold', hjust = 0.5) 

#   legend <- cowplot::get_legend(plots[[1]])

#   for (i in 1:length(plots)) { plots[[i]] <- plots[[i]] + ggtitle(titles[[i]], subtitle = NULL) + theme(legend.position = 'none') } 

# 	plt <- cowplot::plot_grid(la_title, ma_title, NULL, 
#                             plots[[1]], 
# 														plots[[5]], 
# 														NULL,
# 														plots[[2]], 
# 														plots[[6]], 
# 														NULL,
# 														plots[[3]], 
# 														plots[[7]], 
# 														NULL,
# 														plots[[4]], 
# 														plots[[8]],
# 														NULL, ncol = 3,   rel_heights = c(0.15, 1, 1, 1, 1)) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

#   outcome_lookup <- c(prevalence = 'Prevalence per 100,000', incidence = 'Incidence per 100,000', 'Diagnoses per 100,000')

#   y.grob <- textGrob(outcome_lookup, gp=gpar(fontface="bold", fontsize=15), rot=90)

#   x.grob <- textGrob("Year", 
#                    gp=gpar(fontface="bold", fontsize=15))

#   plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

#   ggsave(plt, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), paste0(outcome, '_young.png')))

# }

# make_4x2_old_group_by_state_outcomes_plot <- function(outcome) { 

#   require(cowplot)

#   plots <- list()

#   plots[['la_of']] <- plot_state_outcomes_by_sex('Louisiana', 'oldfemales', outcome)
#   plots[['la_om']] <- plot_state_outcomes_by_sex('Louisiana', 'oldmales', outcome)
#   plots[['la_omsm']] <- plot_state_outcomes_by_sex('Louisiana', 'oldmsm', outcome)
#   plots[['la_omsw']] <- plot_state_outcomes_by_sex('Louisiana', 'oldmsw', outcome)

#   plots[['ma_of']] <- plot_state_outcomes_by_sex('Massachusetts', 'oldfemales', outcome)
#   plots[['ma_om']] <- plot_state_outcomes_by_sex('Massachusetts', 'oldmales', outcome)
#   plots[['ma_omsm']] <- plot_state_outcomes_by_sex('Massachusetts', 'oldmsm', outcome)
#   plots[['ma_omsw']] <- plot_state_outcomes_by_sex('Massachusetts', 'oldmsw', outcome)

# 	la_title <- ggdraw() + 
#     draw_label(
# 			"Louisiana",
# 			fontface = 'bold', hjust = 0.5)

# 	ma_title <- ggdraw() + 
#     draw_label(
# 			"Massachusetts",
# 			fontface = 'bold', hjust = 0.5) 

#   legend <- cowplot::get_legend(plots[[1]])

#   titles <- c('Old F', 
# 							'Old M', 
# 							'Old MSM', 
# 							'Old MSW', 
# 							'Old F', 
# 							'Old M', 
# 							'Old MSM', 
# 							'Old MSW')


#   for (i in 1:length(plots)) { plots[[i]] <- plots[[i]] + ggtitle(titles[[i]], subtitle = NULL) + theme(legend.position = 'none') } 

# 	plt <- cowplot::plot_grid(la_title, ma_title, NULL, 
#                             plots[[1]], 
# 														plots[[5]], 
# 														NULL,
# 														plots[[2]], 
# 														plots[[6]], 
# 														NULL,
# 														plots[[3]], 
# 														plots[[7]], 
# 														NULL,
# 														plots[[4]], 
# 														plots[[8]],
# 														NULL, ncol = 3,   rel_heights = c(0.15, 1, 1, 1, 1)) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

#   outcome_lookup <- c(prevalence = 'Prevalence per 100,000', incidence = 'Incidence per 100,000', 'Diagnoses per 100,000')

#   y.grob <- textGrob(outcome_lookup, gp=gpar(fontface="bold", fontsize=15), rot=90)

#   x.grob <- textGrob("Year", 
#                    gp=gpar(fontface="bold", fontsize=15))

#   plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

#   ggsave(plt, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), paste0(outcome, '_old.png')))

# }


make_incidence_diagnoses_intervention_plots <- function() { 
  # df <- simulate_outcomes_by_sex_both_states()
  df <- readRDS(system.file("interventions/intvs_formatted.rds", package='syphilis'))

  la_inc <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence')
  ma_inc <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'incidence')
  la_diag <- plot_state_outcomes_by_sex_and_age('Louisiana', 'diagnosed')
  ma_diag <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'diagnosed')
  la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence')
  ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence')

  la_inc_w_legend <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_legend=TRUE)
  intervention_legend <- cowplot::get_legend(la_inc_w_legend)

	cowplot::plot_grid(la_diag, ma_diag,  
		la_inc, ma_inc, 
    la_prev, ma_prev,
		nrow=3,  
		labels = c("A", "B", "C", "D", "E", "F")) -> fig4_wout_legend

  cowplot::plot_grid(fig4_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> fig4_new

	ggsave(plot=fig4_new, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'incidence_diagnosis_and_prevalence.png'), units="in", width=12, height=12)

  
  la_inc <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', scales='free_y')
  ma_inc <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'incidence', scales='free_y')
  la_diag <- plot_state_outcomes_by_sex_and_age('Louisiana', 'diagnosed', scales='free_y')
  ma_diag <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'diagnosed', scales='free_y')
  la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence', scales='free_y')
  ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence', scales='free_y')

  la_inc_w_legend <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_legend=TRUE)
  intervention_legend <- cowplot::get_legend(la_inc_w_legend)

	cowplot::plot_grid(la_diag, ma_diag,  
		la_inc, ma_inc, 
    la_prev, ma_prev,
		nrow=3,  
		labels = c("A", "B", "C", "D", "E", "F")) -> fig4_wout_legend

  cowplot::plot_grid(fig4_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> fig4_new

	ggsave(plot=fig4_new, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'incidence_diagnosis_and_prevalence_free_scales.png'), units="in", width=12, height=12)
  # la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence')
  # ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence')

	# cowplot::plot_grid(
		# la_prev, ma_prev, 
		# nrow=1,  
		# labels = c("A", "B")) -> prev_plt_wout_legend

  # cowplot::plot_grid(prev_plt_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> prev_plt

	# ggsave(plot=prev_plt, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'prevalence.png'), units="in", width=12, height=8)

  df %<>% filter(intervention %in% c('basecase', 'msm_annual_hr_msm_quarterly'))

  la_inc <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_ci = T)
  ma_inc <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'incidence', include_ci = T)
  la_diag <- plot_state_outcomes_by_sex_and_age('Louisiana', 'diagnosed', include_ci = T)
  ma_diag <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'diagnosed', include_ci = T)
  la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence', include_ci = T)
  ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence', include_ci = T)

  la_inc_w_legend <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_legend=TRUE)
  intervention_legend <- cowplot::get_legend(la_inc_w_legend)

	cowplot::plot_grid(la_diag, ma_diag,  
		la_inc, ma_inc, 
    la_prev, ma_prev,
		nrow=3,  
		labels = c("A", "B", "C", "D", "E", "F")) -> fig4_wout_legend

  cowplot::plot_grid(fig4_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> fig4_new

	ggsave(plot=fig4_new, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'incidence_diagnosis_and_prevalence_guidelines.png'), units="in", width=12, height=12)

  df %<>% filter(intervention %in% c('basecase', 'msm_annual_hr_msm_quarterly'))

  la_inc <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_ci = T, scales='free_y')
  ma_inc <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'incidence', include_ci = T, scales='free_y')
  la_diag <- plot_state_outcomes_by_sex_and_age('Louisiana', 'diagnosed', include_ci = T, scales='free_y')
  ma_diag <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'diagnosed', include_ci = T, scales='free_y')
  la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence', include_ci = T, scales='free_y')
  ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence', include_ci = T, scales='free_y')

  la_inc_w_legend <- plot_state_outcomes_by_sex_and_age('Louisiana', 'incidence', include_legend=TRUE)
  intervention_legend <- cowplot::get_legend(la_inc_w_legend)

	cowplot::plot_grid(la_diag, ma_diag,  
		la_inc, ma_inc, 
    la_prev, ma_prev,
		nrow=3,  
		labels = c("A", "B", "C", "D", "E", "F")) -> fig4_wout_legend

  cowplot::plot_grid(fig4_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> fig4_new

	ggsave(plot=fig4_new, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'incidence_diagnosis_and_prevalence_guidelines_free_scales.png'), units="in", width=12, height=12)

  
  # la_prev <- plot_state_outcomes_by_sex_and_age('Louisiana', 'prevalence', include_ci = T)
  # ma_prev <- plot_state_outcomes_by_sex_and_age('Massachusetts', 'prevalence', include_ci = T)

	# cowplot::plot_grid(
		# la_prev, ma_prev, 
		# nrow=1,  
		# labels = c("A", "B")) -> prev_plt_wout_legend

  # cowplot::plot_grid(prev_plt_wout_legend, intervention_legend, nrow = 1, rel_widths = c(1,.3)) -> prev_plt

	# ggsave(plot=prev_plt, filename = file.path(system.file('figures/intervention_figures/outcomes_by_age_and_sex/', package = 'syphilis'), 'prevalence_guidelines.png'), units="in", width=12, height=8)
}


save_outcome_plots_with_confidence_intervals <- function() { 
  df <- readRDS(system.file('interventions/intvs_formatted.rds', package='syphilis'))

  options_grid <- expand.grid(
    c('Louisiana', 'Massachusetts'),
    c('prevalence', 'diagnosed', 'incidence'),
    c('youngfemales', 'youngmales', 'oldfemales', 'oldmales'))

  for (iter in 1:nrow(options_grid)) { 
    plt <- plot_state_outcomes_by_sex_and_age_conf_ints(state = options_grid[[iter,1]], outcome = as.character(options_grid[[iter,2]]), group = options_grid[[iter,3]])
    ggsave(plot = plt, filename =
      file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
      package='syphilis'),
      paste0(paste0(sapply(options_grid[iter,], as.character), collapse = '_'), '.png')),
      units = 'in',
      width = 12, height = 8)
  }


  la_prev <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'prevalence', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'prevalence', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'prevalence', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'prevalence', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ma_prev <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'prevalence', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'prevalence', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'prevalence', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'prevalence', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ggsave(plot = la_prev, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'la_prev.png'),
    units = 'in',
    width = 12, height = 8)

  ggsave(plot = ma_prev, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'ma_prev.png'),
    units = 'in',
    width = 12, height = 8)



  la_diag <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'diagnosed', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'diagnosed', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'diagnosed', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'diagnosed', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ma_diag <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'diagnosed', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'diagnosed', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'diagnosed', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'diagnosed', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ggsave(plot = la_diag, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'la_diag.png'),
    units = 'in',
    width = 12, height = 8)

  ggsave(plot = ma_diag, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'ma_diag.png'),
    units = 'in',
    width = 12, height = 8)


  la_inc <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'incidence', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'incidence', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'incidence', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Louisiana', 'incidence', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ma_inc <- cowplot::plot_grid(
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'incidence', group = 'youngfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'incidence', group = 'oldfemales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'incidence', group = 'youngmales'),
    plot_state_outcomes_by_sex_and_age_conf_ints('Massachusetts', 'incidence', group = 'oldmales'),
    nrow = 2, 
    labels = c("A", "B", "C", "D"))

  ggsave(plot = la_inc, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'la_inc.png'),
    units = 'in',
    width = 12, height = 8)

  ggsave(plot = ma_inc, filename =
    file.path(system.file('figures/intervention_figures/with_confidence_intervals/',
    package='syphilis'), 'ma_inc.png'),
    units = 'in',
    width = 12, height = 8)

    
}


plot_point_estimates_of_interventions <- function(state, outcome) { 

  if (! exists('df') || 
  colnames(df) != c(
  "state","intervention","year","variable","group","ci_high","ci_low","mean")) { 
    stop("df should be read in from inst/interventions/intvs_formatted.rds")
  }

  
  df %<>% filter(year == max(year))

  df %<>% mutate(group = recode(group, youngmales = "M: 25-44 y", youngfemales = "F: 25-44 y", oldmales = "M: 45-64 y", oldfemales = "F: 45-64 y"))

  df %<>% mutate(intervention = recode(intervention, 
    basecase = 'Base Case', msm_annual_hr_msm_quarterly = 'Guidelines in MSM',
    msm_quarterly = 'MSM every 3 months', prior_diagnosis_quarterly = 'Prior Diagnosis every 3 months',
    high_activity_quarterly = 'High Activity every 3 months'))


  outcome1 <- enquo(outcome)
  state1 <- enquo(state)

  outcome_lookup <- c(prevalence = "\nPrevalent early syphilis infections per 100,000 population\n", 
    incidence = "\nIncident syphilis cases per 100,000 population\n",
    diagnosed = "\nReported early syphilis cases per 100,000 population\n")

  ggplot(data = filter(df, variable == !! outcome1, state == !! state1), 
    mapping = aes(x = intervention, y = mean, ymax = ci_high, ymin = ci_low, 
    color = intervention, shape = intervention, fill = intervention)) + 
    geom_pointrange() + 
    facet_wrap(~ group) +
    theme_minimal() + 
    coord_flip() + 
    scale_x_discrete(limits = rev(levels(df$intervention))) + 
    xlab("") + 
    ylab(outcome_lookup[[outcome]]) + 
    ggtitle(state) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = .5)) + 
    scale_color_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_fill_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_shape_manual(
      values = c(16,15,17,23,25))


}

plot_point_estimates_panel <- function() { 
  la_diag <- plot_point_estimates_of_interventions('Louisiana', 'diagnosed')
  ma_diag <- plot_point_estimates_of_interventions('Massachusetts', 'diagnosed')
  la_inc <- plot_point_estimates_of_interventions('Louisiana', 'incidence')
  ma_inc <- plot_point_estimates_of_interventions('Massachusetts', 'incidence')
  la_prev <- plot_point_estimates_of_interventions('Louisiana', 'prevalence')
  ma_prev <- plot_point_estimates_of_interventions('Massachusetts', 'prevalence')

  cowplot::plot_grid(
    la_diag, ma_diag, la_inc, ma_inc, nrow = 2, labels = c("A", "B", "C", "D")) -> 
  diagnosis_and_incidence_2016

	ggsave(plot = diagnosis_and_incidence_2016, filename =
	file.path(system.file("figures/intervention_figures/with_confidence_intervals/",
		package='syphilis'),
		'diagnosis_and_incidence_2016.png'),
    units = 'in',
    width = 12, 
    height = 8)

}



plot_point_estimates_of_change_in_interventions <- function(state, outcome) { 

  if (! exists('df') || 
  ! identical(colnames(df),
      c("state", "year", "intervention", "outcome", "ci_high",
      "ci_low", "iq_high", "iq_low", "mean", "median"))
  ) { 
    stop("df should be read in from inst/interventions/intvs_change_from_basecase_summarized.rds")
  }


  # df %<>% mutate(group = recode(group, young_males = "M: 25-44 y", young_females = "F: 25-44 y", 
  #   old_males = "M: 45-64 y", old_females = "F: 45-64 y"))

  df %<>% mutate(intervention = recode(intervention, 
    basecase = 'Base Case', msm_annual_hr_msm_quarterly = 'Guidelines in MSM',
    msm_quarterly = 'MSM every 3 months', prior_diagnosis_quarterly = 'Prior Diagnosis every 3 months',
    high_activity_quarterly = 'High Activity every 3 months'))

  df$intervention %<>% droplevels

  outcome1 <- enquo(outcome)
  state1 <- enquo(state)

  outcome_lookup <- c(prevalence = "\nPercent change in prevalent early syphilis infections\nrelative to the base case\n", 
    incidence = "\nPercent change in incident syphilis cases\nrelative to the base case\n",
    diagnosed = "\nPercent change in reported early syphilis cases\nrelative to the base case\n")

  ggplot(data = filter(df, outcome == !! outcome1, state == !! state1), 
    mapping = aes(x = intervention, middle = median, ymax = ci_high, ymin = ci_low, 
    upper = iq_high, lower = iq_low, 
    color = intervention )) + 
    geom_hline(yintercept = 0, size = .5, color = 'red', alpha = 0.4) +
    # geom_pointrange() + 
    geom_boxplot(stat = 'identity', width = 0.5, alpha = 0.7) + 
    # geom_vline(aes(xintercept = mean)) + 
    # facet_wrap(~ group) +
    theme_minimal() + 
    coord_flip() + 
    expand_limits(y=c(-100, max(df$ci_high, na.rm=T))) + 
    scale_x_discrete(limits = rev(levels(df$intervention))) + 
    xlab("") + 
    ylab(outcome_lookup[[outcome]]) + 
    ggtitle(state) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = .5)) + 
    scale_color_manual(
      labels = c(
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c( "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_fill_manual(
      labels = c( 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c( "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_shape_manual(
      values = c(15,17,23,25))

}


plot_point_estimates_of_relative_change_panel <- function() { 

  df <- readRDS(system.file("interventions/intvs_change_from_basecase_summarized.rds", package='syphilis'))

  la_diag <- plot_point_estimates_of_change_in_interventions('Louisiana', 'diagnosed')
  ma_diag <- plot_point_estimates_of_change_in_interventions('Massachusetts', 'diagnosed')
  la_inc <- plot_point_estimates_of_change_in_interventions('Louisiana', 'incidence')
  ma_inc <- plot_point_estimates_of_change_in_interventions('Massachusetts', 'incidence')
  la_prev <- plot_point_estimates_of_change_in_interventions('Louisiana', 'prevalence')
  ma_prev <- plot_point_estimates_of_change_in_interventions('Massachusetts', 'prevalence')

  cowplot::plot_grid(
    la_diag, ma_diag, la_inc, ma_inc, la_prev, ma_prev, nrow = 3, labels = c("A", "B", "C", "D", "E", "F")) -> 
  inc_diag_and_prev_2017

	ggsave(plot = inc_diag_and_prev_2017, filename =
	file.path(system.file("figures/intervention_figures/with_confidence_intervals/",
		package='syphilis'),
		'change_in_inc_diag_prev_start_of_2017.png'),
    units = 'in',
    width = 12, 
    height = 12)


  # cowplot::plot_grid(
  #   la_prev, ma_prev, nrow = 1, labels = c("A", "B")) -> 
  # prevalence_2017

	# ggsave(plot = prevalence_2017, filename =
	# file.path(system.file("figures/intervention_figures/with_confidence_intervals/",
		# package='syphilis'),
		# 'change_in_prevalence_start_of_2017.png'),
  #   units = 'in',
  #   width = 12, 
  #   height = 8)

}


plot_violin_boxplot_for_change_in_intvs <- function(state, outcome) { 

  if ( ! identical( colnames(df), 
    c("state", "sim", "intervention", "variable", "value"))) { 
      stop("df should be loaded in from prep_data_for_violin_boxplot()") } 

  df %<>% mutate(intervention = recode(intervention, 
    basecase = 'Base Case', msm_annual_hr_msm_quarterly = 'Guidelines in MSM',
    msm_quarterly = 'MSM every 3 months', prior_diagnosis_quarterly = 'Prior Diagnosis every 3 months',
    high_activity_quarterly = 'High Activity every 3 months'))


  df$intervention %<>% droplevels

  outcome1 <- enquo(outcome)
  state1 <- enquo(state)

  outcome_lookup <- c(prevalence = "\nPrevalent early syphilis infections per 100,000 population\n", 
    incidence = "\nIncident syphilis cases per 100,000 population\n",
    diagnosed = "\nReported early syphilis cases per 100,000 population\n")

  ggplot(data = filter(df, variable == !! outcome1, state == !! state1), 
    mapping = aes(x = intervention, y = value, fill = intervention)) + 
    geom_violin(scale = 'width', alpha = 0.8) + 
    geom_boxplot(width = 0.1, color = 'black', fill = 'white', outlier.color = NA) + 
    stat_summary(fun.y=mean, geom='point') + 
    geom_hline(yintercept = 0, size = .5, color = 'red', alpha = 0.4) +
    theme_minimal() + 
    coord_flip() + 
    scale_x_discrete(limits = rev(levels(df$intervention))) + 
    xlab("") + 
    ylab(outcome_lookup[[outcome]]) + 
    ggtitle(state) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = .5)) + 
    scale_color_manual(
      labels = c( 
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) + 
    scale_fill_manual(
      labels = c(
        'Guidelines in MSM', 
        'MSM every 3 months',
        'Prior Diagnosis every 3 months', 
        'High Activity every 3 months'), 
        values = c( "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))

}






plot_violin_boxplot_change_in_intvs_panels <- function() { 

  df <- prep_data_for_violin_boxplot()

  la_diag <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'diagnosed')
  ma_diag <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'diagnosed')
  la_inc <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'incidence')
  ma_inc <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'incidence')
  la_prev <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'prevalence')
  ma_prev <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'prevalence')

  cowplot::plot_grid(
    la_diag, ma_diag, la_inc, ma_inc, la_prev, ma_prev, nrow = 3, labels = c("A", "B", "C", "D", "E", "F")) -> 
  inc_diag_and_prev_2017

	ggsave(plot = inc_diag_and_prev_2017, filename =
	file.path(system.file("figures/intervention_figures/with_confidence_intervals/",
		package='syphilis'),
		'change_in_inc_diag_prev_start_of_2017_violin_version.png'),
    units = 'in',
    width = 12, 
    height = 12)


  df <- prep_data_for_violin_boxplot(exclude_bad_sims = T)


  la_diag <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'diagnosed')
  ma_diag <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'diagnosed')
  la_inc <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'incidence')
  ma_inc <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'incidence')
  la_prev <- plot_violin_boxplot_for_change_in_intvs('Louisiana', 'prevalence')
  ma_prev <- plot_violin_boxplot_for_change_in_intvs('Massachusetts', 'prevalence')

  cowplot::plot_grid(
    la_diag, ma_diag, la_inc, ma_inc, la_prev, ma_prev, nrow = 3, labels = c("A", "B", "C", "D", "E", "F")) -> 
  inc_diag_and_prev_2017

	ggsave(plot = inc_diag_and_prev_2017, filename =
	file.path(system.file("figures/intervention_figures/with_confidence_intervals/",
		package='syphilis'),
		'change_in_inc_diag_prev_start_of_2017_violin_version_exclude_bad_sims.png'),
    units = 'in',
    width = 12, 
    height = 12)
}
