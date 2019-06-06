#' Load Global Variables

#' This function loads the global variables necessary 
#' for the syphLAMA package that do not depend on the 
#' state (LA or MA) having been declared. 

load_globals <- function(model.end = 110) {

	statenames <<- c(LA = 'Louisiana', MA = 'Massachusetts')

  ### set up different subpopulations, sexes, and activity classes ###
  i<<-5 #number of subpopulations, 1= Black, 2=White, 3=Hispanic, 4=MSM-HIVneg, 5=MSM-HIVpos
  j<<-2 # number of sexes, 1=male, 2=female
  k<<-2 #activity class, 1=low, 2=high
  l<<-2 #age groups, 1= 20-44; 2=45-64
  index <<- i*j*k*l

  #### model time steps and calibration period ####
  tstep <<- 1/(52)  #weekly time step
	cal.period <<- 5 #duration of calibration period (2012-2016 currently)
	cal.start <<- 100 #time at which start calibration 
	model.end <<- model.end # end of the simulation
  start.year <<- 2012 # data calibration start
  end.year <<- 2016 # data calibration end
	intervention_years <<- 105:109

	# the term such that cal.start + model_to_gregorian_difference == start.year
	model_to_gregorian_difference <<- start.year - cal.start 

  # background antibiotic treatment
  p.abx.init <- 0.15 #set this to 0 if want to turn off abx treatment -- this implements a period (pre-cal) with high rates of abx use, to represent the intro of penicillin in the pop, which decreases to abx.background 
  p.abx.background <- 0.01
  abx.start <- cal.start-60
  abx.end <- cal.start-20
  p.abx <<- c(rep(0,abx.start), rep(p.abx.init,(abx.end-abx.start)),rep(p.abx.background, (model.end+1-abx.end)) )


  # when to start moving people into the prior infected compartments 
  rep.start <- cal.start - 5 #years before calibration start when start tracking prior treated infections
  rep.count <<- c(rep(0,rep.start), rep(1, model.end+1-rep.start))
  
  age.cat<<-c(25, 20) # age band widths, corresponding to 20-44 yo and 45-64 yo
  out_all <<- NULL
  omega <<- 0.5 #supply and demand of sexual partnerships, 0.5=both sexes copmromise equally, 1=females determine number of partnerhips, 0=males determine number of partnerships

  ### initial model conditions ###
  ind<<- as.matrix(expand.grid(1:k,1:l,1:j,1:i)) #indexing matrix (k,l,j,i)
  init.Y <<-c(rep(c(0,3,0,3),4),rep(c(0,0,0,0),1),rep(c(0,3,0,3),3),rep(c(0,0,0,0),2)) #vector of initial number of infecteds

  #indexes for pulling out different states/population groups from the model output matrix
  s.index <<- 1:index                      #susceptible
  e.index <<- index+1:index                #incubating
  prim.index <<- index*2+1:index           #primary infection
  sec.index <<- index*3+1:index            #secondary infection
  early.index <<- index*4+1:index          #early latent infection
  latent.index <<- index*5+1:index         #late latent infection
  treated.inf.index <<- index*6+1:index    #treated P&S
  treated.early.index <<- index*7+1:index  #treated early latent
  treated.late.index <<- index*8+1:index   #treated late latent
  sr.index <<- index*9+1:index             #re-susceptible
  er.index <<- index*10+1:index            #incubating, previously treated
  primr.index <<- index*11+1:index         #primary infection, previously treated
  secr.index <<- index*12+1:index          #secondary infection, previously treated
  earlyr.index <<- index*13+1:index        #early latent infection, previously treated
  latentr.index <<- index*14+1:index       #late latent infection, previously treated
  inc.index <<- index*15+1:index           #incidence, all
  incr.index <<- index*16+1:index          #incidence, previously treated
  d1.index <<- index*17+1:index            #diagnosed primary
  d2.index <<- index*18+1:index            #diagnosed secondary
  d3.index <<- index*19+1:index            #diagnosed early latent
  nsa.index <<- index*20+1:index           #not sexually active population
  d4.index <<- index*21+1:index            #diagnosed latent latent
  dr.index <<- index*22+1:index            #diagnosed primary, secondary, and early latent, previously treated
	tested.index <<- index*23+1:index        #tested for syphilis

	# index for all infected individuals
	infected.index <<- c(e.index, prim.index, sec.index, early.index,
	  latent.index, er.index, primr.index, secr.index, earlyr.index, latentr.index)
	
	# index for all sexually active individuals
	allpop.index <<- c(infected.index, s.index, treated.inf.index, treated.early.index, treated.late.index, sr.index)

	infectious_index <<- c(prim.index, sec.index, primr.index, secr.index)
	noninfectious_index <<- c(e.index, early.index, latent.index, er.index, earlyr.index, latentr.index)
  
  pop1 <<- c(1:4,21:24) #subpop1
  pop2 <<- c(5:8,25:28) # subpop2
  pop3 <<- c(9:12,29:32) #subpop3
  pop4 <<-c(13:16,33:36 ) #subpop4
  pop5 <<- c(17:20,37:40) #subpop5
  males <<-1:20 # males
  females <<-21:40 #females
  m1<<-1:4 #M subpop1
  m2<<-5:8 #M subpop2
  m3<<-9:12 #M subpop3
  m4<<-13:16 #M subpop4
  m5<<-17:20 #M subpop5
  msw <<-1:12 #M heterosexual
  f1<<-21:24 #F subpop1
  f2<<-25:28 #F subpop2
  f3<<-29:32 #F subpop3
  f4<<-33:36 #F subpop4 (empty)
  f5<<-37:40 #F subpop5 (empty)
  y.m<<-c(1:2,5:6,9:10,13:14,17:18) #youngest age cat M
  y.m.msw<<- c(1:2, 5:6, 9:10) #yougest age cat MSW only
  y.f<<-c(21:22,25:26,29:30) #youngest age cat F
  o.m<<-c(3:4,7:8, 11:12,15:16,19:20) #oldest age cat M
  o.m.msw <<- c(3:4, 7:8, 11:12) #oldest age cat MSW only
  o.f<<-c(23:24,27:28,31:32) #oldest age cat F
  y.m.1 <<- 1:2 #youngest age cat, subpop1
  o.m.1 <<- 3:4 #old age cat, subpop1
  y.m.2 <<- 5:6 #young age cat, subpop 2
  o.m.2 <<- 7:8 # old age cat, subpop 2
  y.m.3 <<- 9:10 # young age cat, subpop 3
  o.m.3 <<- 11:12 # old age cat, subpop 3
  y.m.4 <<- 13:14 # young age cat, subpop 4
  o.m.4 <<- 15:16 # old age cat, subpop 4
  y.m.5 <<- 17:18 #young age cat, subpop 5
  o.m.5 <<- 19:20 #old age cat, subpop 5
  y.f.1 <<- 21:22 #youngest age cat, subpop1
  o.f.1 <<- 23:24 #old age cat, subpop1
  y.f.2 <<- 25:26 #young age cat, subpop 2
  o.f.2 <<- 27:28 # old age cat, subpop 2
  y.f.3 <<- 29:30 # young age cat, subpop 3
  o.f.3 <<- 31:32 # old age cat, subpop 3
  y.f.23 <<- c(25:26,29:30) # youngest age cat, subpops 2&3

  # Make a 4 dimensional array with each of the model dimensions 
	# corresponding to the index of the population in their serial 
	# (1-dimensional) ordering.
  pop_array <<- array(data = 1:40, dim = c(k, l, i, j),
	             dimnames = list(
								k = c('low', 'high'),
								l = c('young', 'old'),
								i = c('black', 'white', 'hispanic', 'msm-hivneg', 'msm-hivpos'),
								j = c('male', 'female')
						   ))

  # Flatten the pop array into a lookup table (data frame)
	# with columns j, k, l, i, index.
	pop <<- as.data.frame.table(pop_array, responseName = 'index')

	# Define screening interventions
	interventions <<- tibble::tribble(
									~codename,                ~longname,      ~target,  ~freq,
										'annual',                 'Annual',        'all',      1,
							'twice_annual',           'Twice Annual',        'all',      2,
								'msm_annual',             'MSM Annual',        'msm',      1,
					'msm_twice_annual',       'MSM Twice Annual',        'msm',      2,
				 'msm_hivpos_annual',        'HIV+ MSM Annual', 'msm-hivpos',      1,
	 'msm_hivpos_twice_annual',  'HIV+ MSM Twice Annual', 'msm-hivpos',      2,
				 'msm_hivneg_annual',        'HIV- MSM Annual', 'msm-hivneg',      1,
	 'msm_hivneg_twice_annual',  'HIV- MSM Twice Annual', 'msm-hivneg',      2
	)
	
	# Natural History Parameters used in Simultaneous Calibration
	natural_history_parameters <<- 
		c('logit.b.m', 'logit.b.f', 'logit.b.msm', 'log.dur.incub', 'log.dur.prim',
			'log.dur.sec', 'log.dur.imm.inf', 'log.dur.imm.early', 'log.dur.immune')

	return(invisible(NULL))
}

load_globals()
