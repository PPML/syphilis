#' Modify the Screening Matrix to Implement Targeted Testing and Treatment Interventions
#' 
#' Check that the intervention is defined appropriately (either not at all or
#' as basecase, or as one of the interventions in the interventions$codename
#' column). Then use a switch on the intervention to call the
#' adjust_screening_for_intervention function with the right arguments.
#' 
modify_simulation_environment_for_an_intervention <- function(e, intervention) {

	if (intervention == 'basecase') return(e)

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
      e$params$screen_repeat[106:nrow(e$params$screen_repeat), ] <- 4
     },

   msm_hr_and_prior_diagnosis_quarterly = {
     adjust_screening_for_intervention(e, pop_indices = c(m4, m5), new_level = 4)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4)
      e$params$screen_repeat[106:nrow(e$params$screen_repeat), ] <- 4
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
      e$params$screen_repeat[106:nrow(e$params$screen_repeat), ] <- 
        4
     },
   high_activity_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 1),
   high_activity_twice_annual = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 2),
   high_activity_quarterly = adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4),
   annual_hr_quarterly = { 
     adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 1)
     adjust_screening_for_intervention(e, pop_indices = high_activity, new_level = 4)
   },
   half_annual = adjust_screening_for_intervention(e, pop_indices = 1:32, new_level = 0.5)
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
adjust_screening_for_intervention <- function(e, intervention_start = 106, pop_indices, new_level, repeat_level = NA) {

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
extractInterventionStatistics <- function(sol) {
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

	data.frame(
		year = intervention_period/52 + model_to_gregorian_difference,
		prevalent_infections = rowSums(sol[1+intervention_period, 1+infected.index]),
		incidents_per_year = rowSums(sol[1+intervention_period, 1+c(inc.index, incr.index)]) - 
			rowSums(sol[1 + intervention_period - 52, 1+c(inc.index, incr.index)]),
		tests_per_year = rowSums(sol[1+intervention_period,1+tested.index]) - 
			rowSums(sol[intervention_period,1+tested.index]),
		popsize = rowSums(sol[1+intervention_period, 1+allpop.index])
		)
}


extract_outcomes_by_sex <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

  decumulate <- function(idxs) { 
    return(
      rowSums(sol[intervention_period, idxs]) - 
      rowSums(sol[intervention_period - 52, idxs])
    )} 

  retrieve <- function(idxs) { 
    rowSums(sol[intervention_period, idxs])
  }

  data.frame(
    year = intervention_period/52 + model_to_gregorian_difference,
    # popsizes
    popsize_all = retrieve(allpop.index),
    popsize_females = retrieve(females),
    popsize_msw = retrieve(msw),
    popsize_msm = retrieve(msm),

    # diagnoses 
    diagnosed_all = decumulate(c(d1.index, d2.index, d3.index)), # diagnosed early syphilis
    diagnosed_females = decumulate(diagnosed_females),
    diagnosed_msw = decumulate(diagnosed_msw),
    diagnosed_msm = decumulate(diagnosed_msm),

    # incidence
    incidence_all = decumulate(c(inc.index, incr.index)),
    incidence_females = decumulate(incidence_females),
    incidence_msw = decumulate(incidence_msw),
    incidence_msm = decumulate(incidence_msm),

    # prevalence
    prevalence_all = retrieve(early_infected_index),
    prevalence_females = retrieve(prevalence_early_females),
    prevalence_msw = retrieve(prevalence_early_msw),
    prevalence_msm = retrieve(prevalence_early_msm),

    # tests
    tests_all = decumulate(tested.index)
  )

}


extract_diagnoses_by_sex <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

  decumulate <- function(idxs) { 
    return(
      rowSums(sol[intervention_period, idxs]) - 
      rowSums(sol[intervention_period - 52, idxs])
    )} 

  data.frame(
    year = intervention_period/52 + model_to_gregorian_difference,
    diagnosed_females = decumulate(diagnosed_females),
    diagnosed_msw = decumulate(diagnosed_msw),
    diagnosed_msm = decumulate(diagnosed_msm),
    popsize_females = rowSums(sol[intervention_period, females]),
    popsize_msw = rowSums(sol[intervention_period, msw]),
    popsize_msm = rowSums(sol[intervention_period, msm])
  )
}


extract_incidence_by_sex <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

  decumulate <- function(idxs) { 
    return(
      rowSums(sol[intervention_period, idxs]) - 
      rowSums(sol[intervention_period - 52, idxs])
    )} 

  data.frame(
    year = intervention_period/52 + model_to_gregorian_difference,
    incidence_females = decumulate(incidence_females),
    incidence_msw = decumulate(incidence_msw),
    incidence_msm = decumulate(incidence_msm),
    popsize_females = rowSums(sol[intervention_period, females]),
    popsize_msw = rowSums(sol[intervention_period, msw]),
    popsize_msm = rowSums(sol[intervention_period, msm])
  )
}


extract_prevalence_by_sex <- function(sol) { 
	intervention_start <- min(intervention_years)
	intervention_period <- seq((intervention_start-1)*52-1, model.end*52)

  decumulate <- function(idxs) { 
    return(
      rowSums(sol[intervention_period, idxs]) - 
      rowSums(sol[intervention_period - 52, idxs])
    )} 

  data.frame(
    year = intervention_period/52 + model_to_gregorian_difference,
    prevalence_females = rowSums(sol[intervention_period, prevalence_females]),
    prevalence_msw = rowSums(sol[intervention_period, prevalence_msw]),
    prevalence_msm = rowSums(sol[intervention_period, prevalence_msm]),
    popsize_females = rowSums(sol[intervention_period, females]),
    popsize_msw = rowSums(sol[intervention_period, msw]),
    popsize_msm = rowSums(sol[intervention_period, msm])
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
simulate_interventions <- function(
  theta, 
  interventions_list = c('msm_annual_hr_msm_quarterly',
  'prior_diagnosis_quarterly', 'high_activity_quarterly',
  'msm_annual_hr_and_prior_diagnosis_quarterly',
  'annual',
  'half_annual',
  'msm_hr_and_prior_diagnosis_quarterly',
  'msm_annual_all_hr_half_annual',
  'msm_annual_all_hr_annual',
  'twelve_times_annual',
  'sixteen_times_annual',
  'twenty_four_times_annual'
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
  df <- simulate_interventions(theta, info_extractor = extract_outcomes_by_sex)

  df %<>% 
  # reshape data into rates per 100k
  dplyr::mutate(
    # prevalence 
    prevalence_all = prevalence_all / popsize_all * 1e5,
    prevalence_females = prevalence_females / popsize_females * 1e5,
    prevalence_msw = prevalence_msw / popsize_msw * 1e5,
    prevalence_msm = prevalence_msm / popsize_msm * 1e5,

    # incidence 
    incidence_all = incidence_all / popsize_all * 1e5,
    incidence_females = incidence_females / popsize_females * 1e5,
    incidence_msw = incidence_msw / popsize_msw * 1e5,
    incidence_msm = incidence_msm / popsize_msm * 1e5,

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
  tidyr::separate(col = variable, into = c('variable', 'sex'), extra='drop')

  return(df)
}


simulate_diagnoses_by_sex_by_state <- function(state) { 
  
  state <<- state
  # load state model
  load.start()

  # get state parametrization
  theta <- get(paste0('theta_', tolower(state)))

  # simulate interventions, record diagnoses by sex
  df <- simulate_interventions(theta, info_extractor = extract_diagnoses_by_sex)

  # reshape data into diagnoses per 100k
  df %<>% 
  mutate(diagnosed_females_per_100k = diagnosed_females / popsize_females * 100000,
         diagnosed_msw_per_100k = diagnosed_msw / popsize_msw * 100000,
         diagnosed_msm_per_100k = diagnosed_msm / popsize_msm * 100000) %>% 
  select(-c(diagnosed_females, diagnosed_msw, diagnosed_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
  reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'diagnosis_rate') %>% 
  tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
  select(-prefix)

  return(df)
}


simulate_incidence_by_sex_by_state <- function(state) { 
  
  state <<- state
  # load state model
  load.start()

  # get state parametrization
  theta <- get(paste0('theta_', tolower(state)))

  # simulate interventions, record diagnoses by sex
  df <- simulate_interventions(theta, info_extractor = extract_incidence_by_sex)

  # reshape data into diagnoses per 100k
  df %<>% 
  mutate(incidence_females_per_100k = incidence_females / popsize_females * 100000,
         incidence_msw_per_100k = incidence_msw / popsize_msw * 100000,
         incidence_msm_per_100k = incidence_msm / popsize_msm * 100000) %>% 
  select(-c(incidence_females, incidence_msw, incidence_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
  reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'incidence_rate') %>% 
  tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
  select(-prefix)

  return(df)
}


simulate_prevalence_by_sex_by_state <- function(state) { 
  
  state <<- state
  # load state model
  load.start()

  # get state parametrization
  theta <- get(paste0('theta_', tolower(state)))

  # simulate interventions, record diagnoses by sex
  df <- simulate_interventions(theta, info_extractor = extract_prevalence_by_sex)

  # reshape data into diagnoses per 100k
  df %<>% 
  mutate(prevalence_females_per_100k = prevalence_females / popsize_females * 100000,
         prevalence_msw_per_100k = prevalence_msw / popsize_msw * 100000,
         prevalence_msm_per_100k = prevalence_msm / popsize_msm * 100000) %>% 
  select(-c(prevalence_females, prevalence_msw, prevalence_msm, popsize_females, popsize_msw, popsize_msm)) %>% 
  reshape2::melt(id.vars = c('intervention', 'year'), value.name = 'prevalence_rate') %>% 
  tidyr::separate(col = variable, into = c('prefix', 'group'), extra='drop') %>% 
  select(-prefix)

  return(df)
}

simulate_outcomes_by_sex_both_states <- function() { 
  rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', simulate_outcomes_by_sex_by_state('LA')),
    cbind.data.frame(state = 'Massachusetts', simulate_outcomes_by_sex_by_state('MA')))
}

#' @example 
#' df <- simulate_outcomes_by_sex_both_states()
#' plot_list <- plot_outcomes_by_sex(df, save=F)
plot_outcomes_by_sex <- function(df) { 
  plot_list <- list()

  df %<>% filter(year <= 2027, year >= 2017)

  produce_plot <- function(state, sex, outcome, include_title = T) { 
    state1 <- enquo(state)
    sex1 <- enquo(sex)
    outcome1 <- enquo(outcome)

    sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM', all = 'All')

    ggplot( data = filter(df, state == !! state1, sex == !! sex1, variable == !! outcome1), 
            mapping = aes(x = year, y = rate, color = intervention, linetype = intervention)) + 
      geom_line(size = 1) + 
      ylab('') + 
      xlab('') + 
      expand_limits(y=0) + 
      # scale_color_manual(values = c('red', 'forestgreen', 'skyblue', 'orange', 'coral', 'purple', 'pink', 'grey', 'grey1', 'grey2', 'grey3', 'grey4')) + 
      scale_x_continuous(breaks = seq(2017, 2027, by=2)) + 
      theme_bw() + 
      if (include_title) { 
        ggtitle(paste0(tools::toTitleCase(outcome), ' Rate per 100,000'), subtitle = paste0(sex_lookup[[sex]], " - ", state)) 
      } else { NULL } 
  }
}



plot_state_outcomes_by_sex <- function(state, sex, outcome, include_title = T) { 

  if (! exists('df') || 
  colnames(df) != c("state", "intervention", "year", "variable", "sex", "rate")) { 
    stop("df should be produced using simulate_outcomes_by_sex_both_states()")
  }

  state1 <- enquo(state)
  sex1 <- enquo(sex)
  outcome1 <- enquo(outcome)

  sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM', all = 'All')

  ggplot( data = filter(df, state == !! state1, sex == !! sex1, variable == !! outcome1), 
          mapping = aes(x = year, y = rate, color = intervention, linetype = intervention)) + 
    geom_line(size = 1) + 
    ylab('') + 
    xlab('') + 
    expand_limits(y=0) + 
    # scale_color_manual(values = c('red', 'forestgreen', 'skyblue', 'orange', 'coral', 'purple', 'pink', 'grey', 'grey1', 'grey2', 'grey3', 'grey4')) + 
    scale_x_continuous(breaks = seq(2017, 2027, by=2)) + 
    theme_bw() + 
    if (include_title) { 
      ggtitle(paste0(tools::toTitleCase(outcome), ' Rate per 100,000'), subtitle = paste0(sex_lookup[[sex]], " - ", state)) 
    } else { NULL } 
}

#' @import tidyverse
simulate_and_plot_diagnoses_by_sex_both_states <- function() { 
	load_globals(model.end = 116)

  la_df <- simulate_diagnoses_by_sex_by_state('LA')
  ma_df <- simulate_diagnoses_by_sex_by_state('MA')

  df <- rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', la_df),
    cbind.data.frame(state = 'Massachusetts', ma_df))

  df %<>% filter(year <= 2027, year >= 2017)

  produce_plot <- function(state, sex) { 
    state1 <- enquo(state)
    sex1 <- enquo(sex)

    sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

    ggplot( data = filter(df, state == !! state, group == !! sex), 
            mapping = aes(x = year, y = diagnosis_rate, color = intervention, linetype = intervention)) + 
      geom_line(size = 1) + 
      ylab('') + 
      xlab('') + 
      expand_limits(y=0) + 
      scale_color_manual(values = c('red', 'forestgreen', 'skyblue', 'orange', 'coral')) + 
      scale_x_continuous(breaks = seq(2017, 2027, by=2)) + 
      theme_bw() + 
      ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
  }


  p1 <- produce_plot('Louisiana', 'females') 
  p2 <- produce_plot('Massachusetts', 'females')
  p3 <- produce_plot('Louisiana', 'msw')
  p4 <- produce_plot('Massachusetts', 'msw')
  p5 <- produce_plot('Louisiana', 'msm')
  p6 <- produce_plot('Massachusetts', 'msm')

  legend <- cowplot::get_legend(p1)

  p1 <- p1 + theme(legend.position = 'none')
  p2 <- p2 + theme(legend.position = 'none')
  p3 <- p3 + theme(legend.position = 'none')
  p4 <- p4 + theme(legend.position = 'none')
  p5 <- p5 + theme(legend.position = 'none')
  p6 <- p6 + theme(legend.position = 'none')

  plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

  y.grob <- textGrob("Diagnosed Cases Per 100,000", 
                   gp=gpar(fontface="bold", fontsize=15), rot=90)

  x.grob <- textGrob("Year", 
                   gp=gpar(fontface="bold", fontsize=15))

  plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

  # ggsave(plot = plt, 
  #   filename = file.path(system.file("figures/diagnoses_by_sex_in_interventions/", package = 'syphilis'), 'diagnoses_by_sex_in_interventions.png'),
  #   width = 10)

  return(plt)
}


#' @import tidyverse
simulate_and_plot_incidence_by_sex_both_states <- function() { 
	load_globals(model.end = 116)

  la_df <- simulate_incidence_by_sex_by_state('LA')
  ma_df <- simulate_incidence_by_sex_by_state('MA')

  df <- rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', la_df),
    cbind.data.frame(state = 'Massachusetts', ma_df))

  df %<>% filter(year <= 2027, year >= 2017)

  produce_plot <- function(state, sex) { 
    state1 <- enquo(state)
    sex1 <- enquo(sex)

    sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

    ggplot( data = filter(df, state == !! state, group == !! sex), 
            mapping = aes(x = year, y = incidence_rate, color = intervention, linetype = intervention)) + 
      geom_line(size = 1) + 
      ylab('') + 
      xlab('') + 
      expand_limits(y=0) + 
      scale_color_manual(
        labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'), 
        values = c('red', 'forestgreen', 'skyblue', 'orange', 'purple')) + 
      scale_linetype_manual(
        labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'),
        values = c(1,5,2,3,4))  +
      scale_x_continuous(breaks = seq(2017, 2027, by=2)) + 
      theme_bw() + 
      
      ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
  }


  p1 <- produce_plot('Louisiana', 'females')
  p2 <- produce_plot('Massachusetts', 'females')
  p3 <- produce_plot('Louisiana', 'msw')
  p4 <- produce_plot('Massachusetts', 'msw')
  p5 <- produce_plot('Louisiana', 'msm')
  p6 <- produce_plot('Massachusetts', 'msm')

  legend <- cowplot::get_legend(p1)

  p1 <- p1 + theme(legend.position = 'none')
  p2 <- p2 + theme(legend.position = 'none')
  p3 <- p3 + theme(legend.position = 'none')
  p4 <- p4 + theme(legend.position = 'none')
  p5 <- p5 + theme(legend.position = 'none')
  p6 <- p6 + theme(legend.position = 'none')

  plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

  y.grob <- textGrob("Incidence Per 100,000", 
                   gp=gpar(fontface="bold", fontsize=15), rot=90)

  x.grob <- textGrob("Year", 
                   gp=gpar(fontface="bold", fontsize=15))

  plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

  ggsave(plot = plt, 
    filename = file.path(system.file("figures/", package = 'syphilis'), 'incidence_by_sex_in_interventions.png'),
    width = 10)

  return(plt)
}

#' @import tidyverse
simulate_and_plot_prevalence_by_sex_both_states <- function() { 
	load_globals(model.end = 116)

  la_df <- simulate_prevalence_by_sex_by_state('LA')
  ma_df <- simulate_prevalence_by_sex_by_state('MA')

  df <- rbind.data.frame(
    cbind.data.frame(state = 'Louisiana', la_df),
    cbind.data.frame(state = 'Massachusetts', ma_df))

  df %<>% filter(year <= 2027, year >= 2017)

  produce_plot <- function(state, sex) { 
    state1 <- enquo(state)
    sex1 <- enquo(sex)

    sex_lookup <- c(females = 'F', msw = 'MSW', msm = 'MSM')

    ggplot( data = filter(df, state == !! state, group == !! sex), 
            mapping = aes(x = year, y = prevalence_rate, color = intervention, linetype = intervention)) + 
      geom_line(size = 1) + 
      ylab('') + 
      xlab('') + 
      expand_limits(y=0) + 
      scale_color_manual(
        labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'), 
        values = c('red', 'forestgreen', 'skyblue', 'orange', 'purple')) + 
      scale_linetype_manual(
        labels = c('Base Case', 'Guidelines in MSM', 'Prior Diagnosis Quarterly', 'High Activity Quarterly', 'Combination'),
        values = c(1,5,2,3,4))  +
      scale_x_continuous(breaks = seq(2017, 2027, by=2)) + 
      theme_bw() + 
      
      ggtitle(paste0(sex_lookup[[sex]], " - ", state)) 
  }


  p1 <- produce_plot('Louisiana', 'females')
  p2 <- produce_plot('Massachusetts', 'females')
  p3 <- produce_plot('Louisiana', 'msw')
  p4 <- produce_plot('Massachusetts', 'msw')
  p5 <- produce_plot('Louisiana', 'msm')
  p6 <- produce_plot('Massachusetts', 'msm')

  legend <- cowplot::get_legend(p1)

  p1 <- p1 + theme(legend.position = 'none')
  p2 <- p2 + theme(legend.position = 'none')
  p3 <- p3 + theme(legend.position = 'none')
  p4 <- p4 + theme(legend.position = 'none')
  p5 <- p5 + theme(legend.position = 'none')
  p6 <- p6 + theme(legend.position = 'none')

  plt <- cowplot::plot_grid(p1, p2, NULL, p3, p4, NULL, p5, p6, NULL, ncol = 3) + cowplot::draw_grob(legend, 2/3.3, 0, 1.5/3.3, 1)

  y.grob <- textGrob("Prevalence Per 100,000", 
                   gp=gpar(fontface="bold", fontsize=15), rot=90)

  x.grob <- textGrob("Year", 
                   gp=gpar(fontface="bold", fontsize=15))

  plt <- grid.arrange(arrangeGrob(plt, left = y.grob, bottom = x.grob))

  ggsave(plot = plt, 
    filename = file.path(system.file("figures/", package = 'syphilis'), 'prevalence_by_sex_in_interventions.png'),
    width = 10)

  return(plt)
}

