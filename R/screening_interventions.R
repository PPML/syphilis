#' Modify the Screening Matrix to Implement Targeted Testing and Treatment Interventions
#' 
#' Check that the intervention is defined appropriately (either not at all or
#' as basecase, or as one of the interventions in the interventions$codename
#' column). Then use a switch on the intervention to call the
#' adjust_screening_for_intervention function with the right arguments.
#' 
modify_simulation_environment_for_an_intervention <- function(e, intervention) {

	if (intervention == 'basecase') return(e)

  first_intervention_year <- min(intervention_years)  # +1 for C++ indexing

  # turn the intervention codename into arguments for adjust_screening_for_intervention
	switch(intervention, 
									annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 1),
						twice_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 2),
			 twelve_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 12),
   twenty_four_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 24),
   sixteen_times_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 16),
							msm_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 1),
				msm_twice_annual = adjust_screening_for_intervention(e, pop_indices = c(m4,m5), new_level = 2),
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
	need_increase <- e$screen[intervention_years, pop_indices] < new_level
	e$screen[intervention_years, pop_indices][need_increase] <- new_level
  
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
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

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
    popsize_females = retrieve(all_females),
    popsize_msw = retrieve(all_msw),
    popsize_msm = retrieve(all_msm),

    # diagnoses 
    diagnosed_all = decumulate(c(d1.index, d2.index, d3.index)), # diagnosed early syphilis
    diagnosed_females = decumulate(diagnosed_females),
    diagnosed_msw = decumulate(diagnosed_msw),
    diagnosed_msm = decumulate(diagnosed_msm),
    diagnosed_reinfected = decumulate(dr.index),
    

    # incidence
    incidence_all = decumulate(c(inc.index, incr.index)),
    incidence_females = decumulate(incidence_females),
    incidence_msw = decumulate(incidence_msw),
    incidence_msm = decumulate(incidence_msm),
    incidence_reinfected = decumulate(incr.index),
    incidence_new = decumulate(inc.index),
    incidence_high_activity = decumulate(incidence_high_activity),
    incidence_low_activity = decumulate(incidence_low_activity),



    # prevalence
    prevalence_all = retrieve(early_infected_index),
    prevalence_females = retrieve(prevalence_early_females),
    prevalence_msw = retrieve(prevalence_early_msw),
    prevalence_msm = retrieve(prevalence_early_msm),
    prevalence_reinfected = retrieve(reinfected.index),
    prevalence_new = retrieve(newly_infected.index),
    prevalence_high_activity = retrieve(prevalence_high_activity),
    prevalence_low_activity = retrieve(prevalence_low_activity),
    prevalence_noninfectious = retrieve(noninfectious_index),
    prevalence_infectious = retrieve(infectious_index),

    # tests
    tests_all = decumulate(tested.index)
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
  'prior_diagnosis_annual', 
  'prior_diagnosis_quarterly', 
  'high_activity_annual',
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
    prevalence_all = prevalence_all / popsize_all * 1e5,
    prevalence_females = prevalence_females / popsize_females * 1e5,
    prevalence_msw = prevalence_msw / popsize_msw * 1e5,
    prevalence_msm = prevalence_msm / popsize_msm * 1e5,
    prevalence_new = prevalence_new / popsize_all * 1e5,
    prevalence_reinfected = prevalence_reinfected / popsize_all * 1e5,
    prevalence_high_activity = prevalence_high_activity / popsize_all * 1e5,
    prevalence_low_activity = prevalence_low_activity / popsize_all * 1e5,

    # incidence 
    incidence_all = incidence_all / popsize_all * 1e5,
    incidence_females = incidence_females / popsize_females * 1e5,
    incidence_msw = incidence_msw / popsize_msw * 1e5,
    incidence_msm = incidence_msm / popsize_msm * 1e5,
    incidence_reinfected = incidence_reinfected / popsize_all * 1e5,
    incidence_new = incidence_new / popsize_all * 1e5,
    incidence_high_activity = incidence_high_activity / popsize_all * 1e5,
    incidence_low_activity = incidence_low_activity / popsize_all * 1e5,

    # diagnoses
    diagnosed_all = diagnosed_all / popsize_all * 1e5,
    diagnosed_females = diagnosed_females / popsize_females * 1e5,
    diagnosed_msw = diagnosed_msw / popsize_msw * 1e5,
    diagnosed_msm = diagnosed_msm / popsize_msm * 1e5) %>% 

  # remove population sizes
  select(
    -c(popsize_females, 
       popsize_msw, 
       popsize_msm,
       tests_all)) %>% 

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

plot_state_outcomes_by_sex <- function(state, group, outcome, include_title = T, years) { 

  if (missing(years)) {
    min_year <- intervention_start_gregorian
    max_year <- intervention_stop_gregorian
  }

  if (! exists('df') || 
  colnames(df) != c("state", "intervention", "year", "variable", "group", "rate")) { 
    stop("df should be produced using simulate_outcomes_by_sex_both_states()")
  }

  df %<>% filter(year <= max_year, year >= min_year)

  state1 <- enquo(state)
  group1 <- enquo(group)
  outcome1 <- enquo(outcome)

  group_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM', all = 'All', 
    new = 'New Infections', reinfected = 'Reinfected')

  outcome_lookup <- c(prevalence = "Prevalence of Early Syphilis", incidence = "Incidence Rate", diagnosed = "Diagnosed Cases")

  ggplot(
    data = filter(df, state == !! state1, group == !! group1, variable == !! outcome1), 
    mapping = aes(x = year, y = rate, color = intervention, linetype = intervention)) + 
    geom_line(size = 1.25, alpha=0.8) + 
    ylab('') + 
    xlab('') + 
    expand_limits(y=0) + 
    scale_color_manual(
      labels = c('Base Case', 
        'Guidelines in MSM', 
        'Prior Diagnosis Annual', 
        'Prior Diagnosis Quarterly', 
        'High Activity Annual',
        'High Activity Quarterly'), 
        # values = c("#05141E", "#762B19", "#3D507A", "#657062", "#D14E3E", "#E78A40", "#EBD799")
        # values = c("#7BA46C", "#602D31", "#008D91", "#0A789F", "#C6A28A", "#61B8D3", "#EACF9E")
        # values = c("#7FC97F", "#BEAED4", "#FDC086", "#AAAA99", "#386CB0", "#F0027F", "#BF5B17")
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#AAAA33")) + 

      scale_linetype_manual(
        labels = c('Base Case', 
          'Guidelines in MSM', 
          'Prior Diagnosis Annual', 
          'Prior Diagnosis Quarterly', 
          'High Activity Annual',
          'High Activity Quarterly'), 
        values = c(1,5,2,3,4,6,7))  +
      scale_x_continuous(breaks = seq(intervention_start_gregorian,
        intervention_stop_gregorian, by=2)) + 
      theme_bw() + 
      if (include_title) { 
        ggtitle(paste0(outcome_lookup[[outcome]], ' per 100,000'), subtitle = paste0(group_lookup[[group]], " - ", state)) 
      } else { NULL } 
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


make_intervention_analysis_plots <- function() { 
  df <- simulate_outcomes_by_sex_both_states()

  # df %<>% filter((intervention == 'basecase' & year %% 1 == 0 | year < 2012) | intervention != 'basecase')

  la_prevalence <- plot_state_outcomes_by_sex('Louisiana', 'all', 'prevalence')

  ma_prevalence <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'prevalence')

  ggsave(plot=la_prevalence, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence.png'),
    width = 10, height = 6)


  ggsave(plot=ma_prevalence, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_prevalence.png'),
    width = 10, height = 6)


  la_incidence <- plot_state_outcomes_by_sex('Louisiana', 'all', 'incidence')
  ma_incidence <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'incidence')

  ggsave(plot=la_incidence, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence.png'),
    width = 10, height = 6)

  ggsave(plot=ma_incidence, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence.png'),
    width = 10, height = 6)

  la_diagnosed <- plot_state_outcomes_by_sex('Louisiana', 'all', 'diagnosed')
  ma_diagnosed <- plot_state_outcomes_by_sex('Massachusetts', 'all', 'diagnosed')

  ggsave(plot=la_diagnosed, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_diagnosed.png'),
    width = 10, height = 6)

  ggsave(plot=ma_diagnosed, 
    file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_diagnosed.png'),
    width = 10, height = 6)

  df <- breakdown_basecase()

  plot_breakdown_by_state_outcome_grouping <- function(state, outcome, grouping) { 
    state1 <- enquo(state)
    outcome1 <- enquo(outcome)
    grouping1 <- enquo(grouping)

    grouping_lookup <- list('risk' = c('high', 'low'), sex = c('msm', 'msw', 'females'), reinfection = c('new', 'reinfected'), infectiousness = c('noninfectious', 'infectious'))
    

    grouping_color_lookup <- c(risk = 'Reds', reinfection = 'Purples', sex = 'Blues', infectiousness = 'Oranges')

    ggplot(
        data = filter(df, year >= 2012, year <= 2016, state == !! state1, variable == !! outcome1, group %in% grouping_lookup[[grouping]]),
        mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
        geom_area(alpha=0.8) + 
        scale_fill_brewer(palette = grouping_color_lookup[[grouping]]) + 
        theme_bw() + 
        scale_x_continuous(breaks = seq(2012, 2016, by=2)) + 
        ggtitle(paste0(tools::toTitleCase(outcome), " per 100,000 in ", state), subtitle = paste0("Broken down by ", tools::toTitleCase(grouping)))
  }

  la_prevalence_breakdown_risk <- ggplot(
    data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('high', 'low')),
    mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
    geom_area(alpha=0.5) + 
    scale_fill_brewer(palette = 'Reds') + 
    theme_bw()

  la_prevalence_breakdown_infectiousness <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'prevalence', 'infectiousness')
  ggsave(plot = la_prevalence_breakdown_infectiousness, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence_breakdown_infectiousness.png'))

  ma_prevalence_breakdown_infectiousness <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'prevalence', 'infectiousness')
  ggsave(plot = ma_prevalence_breakdown_infectiousness, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_prevalence_breakdown_infectiousness.png'))

  la_incidence_breakdown_risk <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'risk')
  ggsave(plot = la_incidence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_risk.png'))

  ma_incidence_breakdown_risk <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'risk')
  ggsave(plot = ma_incidence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_risk.png'))

  la_incidence_breakdown_sex <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'sex')
  ggsave(plot = la_incidence_breakdown_sex, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_sex.png'))

  ma_incidence_breakdown_sex <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'sex')
  ggsave(plot = ma_incidence_breakdown_sex, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_sex.png'))


  la_incidence_breakdown_reinfection <- plot_breakdown_by_state_outcome_grouping('Louisiana', 'incidence', 'reinfection')
  ggsave(plot = la_incidence_breakdown_reinfection, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_incidence_breakdown_reinfection.png'))

  ma_incidence_breakdown_reinfection <- plot_breakdown_by_state_outcome_grouping('Massachusetts', 'incidence', 'reinfection')
  ggsave(plot = ma_incidence_breakdown_reinfection, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'ma_incidence_breakdown_reinfection.png'))
  # la_prevalence_breakdown_sex <- ggplot(
  #   data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('msm', 'msw', 'females')),
  #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
  #   geom_area(alpha=0.5) + 
  #   scale_fill_brewer(palette = 'Reds') + 
  #   theme_bw()

  # la_prevalence_breakdown_reinfection <- ggplot(
  #   data = filter(df, state == 'Louisiana', variable == 'prevalence', group %in% c('new', 'reinfected')),
  #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
  #   geom_area(alpha=0.5) + 
  #   scale_fill_brewer(palette = 'Purples') + 
  #   theme_bw()

  # ma_prevalence_breakdown_risk <- ggplot(
  #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('high', 'low')),
  #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
  #   geom_area(alpha=0.5) + 
  #   scale_fill_brewer(palette = 'Reds') + 
  #   theme_bw()

  # ma_prevalence_breakdown_sex <- ggplot(
  #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('msm', 'msw', 'females')),
  #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
  #   geom_area(alpha=0.5) + 
  #   scale_fill_brewer(palette = 'Reds') + 
  #   theme_bw()

  # ma_prevalence_breakdown_reinfection <- ggplot(
  #   data = filter(df, state == 'Massachusetts', variable == 'prevalence', group %in% c('new', 'reinfected')),
  #   mapping = aes(x = year, y = rate, group = group,  fill = group)) + 
  #   geom_area(alpha=0.5) + 
  #   scale_fill_brewer(palette = 'Purples') + 
  #   theme_bw()


  ggsave(plot = la_prevalence_breakdown_risk, filename = file.path(system.file('figures/intervention_figures/', package='syphilis'), 'la_prevalence_breakdown_risk.png'))


}
