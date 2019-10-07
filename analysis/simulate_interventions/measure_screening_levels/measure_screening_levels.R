
devtools::load_all()

measure_screening_levels <- function() {

  outcomes <- list()

  for (state in c('LA', 'MA')) { 
    state <<- state
    load.start()

    theta <- get(paste0('theta_', tolower(state)))

    for (intv in c('basecase', 'msm_annual_hr_msm_quarterly',
      'msm_quarterly', 
      'hr_quarterly', 'prior_diagnosis_quarterly')) {

      e <- constructSimulationEnvironment(theta)

      modify_simulation_environment_for_an_intervention(e, intv)

      runSimulation(e)

      out <- e$sim$out

      total_tests <- sum(out[2017 - model_to_gregorian_difference,tested.index]) -
        sum(out[2012 - model_to_gregorian_difference,tested.index]) 

      outcomes[[length(outcomes)+1]] <- 
        c(state = state, intv = intv, 
        total_tests = total_tests)
    }
  }
  df <- do.call(rbind.data.frame, outcomes)
  colnames(df) <- c('state', 'intv', 'total_tests')
  return(df)
}

df <- measure_screening_levels()

saveRDS(df, file=file.path(system.file('interventions/',
  package='syphilis'), 'intervention_tests.rds'))
write.csv(df, file=file.path(system.file('interventions/',
  package='syphilis'), 'intervention_tests.csv'))

library(dplyr)

basecase <- filter(df, intv == 'basecase')

df <- filter(df, intv != 'basecase')

df2 <- merge(df, basecase, by = c('state'))

df2$total_tests.x %<>% as.character %>% as.numeric
df2$total_tests.y %<>% as.character %>% as.numeric

df2 %<>% mutate(change_in_tests_pct =
  (total_tests.x - total_tests.y)/(total_tests.y) * 100 )

saveRDS(df2, file=file.path(system.file('interventions/',
  package='syphilis'), 'intervention_tests_pct_chng.rds'))

write.csv(df2,file=file.path(system.file('interventions/',
  package='syphilis'), 'intervention_tests_pct_chng.csv'))
 

