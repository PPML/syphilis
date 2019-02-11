#############################################
###   function to calculate aging rates   ###
#############################################
## aging function needs to take into account changing size of non-sexually active population with age
## and maintaining the size of the high and low risk sexual activity groups and HIV+ population

aging.fun <- function(age.cat){
  ages <-rep(c(rep(age.cat[1],2),rep(age.cat[2],2)),j*i) # age band widths, corresponding to 20-44 yo and 45-64 yo, by AC 
  aging <- diag(-1/ages)
  p<-matrix(0,nrow=index, ncol=index) ## rate of entry into older age group is split between different sized ACs
  p[3,1:2]=p.k[1,1,2,1]/age.cat[1] #males
  p[4,1:2]=p.k[1,2,2,1]/age.cat[1]
  p[7,5:6]=p.k[1,1,2,2]/age.cat[1]
  p[8,5:6]=p.k[1,2,2,2]/age.cat[1]
  p[11,9:10]=p.k[1,1,2,3]/age.cat[1]
  p[12,9:10]=p.k[1,2,2,3]/age.cat[1]
  p[15,13:14]=(p.k[1,1,2,4]/age.cat[1])*(1-aging.rate.hiv) #move some younger HIV- MSM into HIV+ group as age, to reflect higher HIV prevalence in older age cat
  p[16,13:14]=(p.k[1,2,2,4]/age.cat[1])*(1-aging.rate.hiv)
  p[19,13:14]=(p.k[1,1,2,5]/age.cat[1])*aging.rate.hiv
  p[20,13:14]=(p.k[1,2,2,5]/age.cat[1])*aging.rate.hiv
  p[19,17:18]=p.k[1,1,2,5]/age.cat[1]
  p[20,17:18]=p.k[1,2,2,5]/age.cat[1]
  
  p[23,21:22]=p.k[2,1,2,1]/age.cat[1] #females
  p[24,21:22]=p.k[2,2,2,1]/age.cat[1]
  p[27,25:26]=p.k[2,1,2,2]/age.cat[1]
  p[28,25:26]=p.k[2,2,2,2]/age.cat[1]
  p[31,29:30]=p.k[2,1,2,3]/age.cat[1]
  p[32,29:30]=p.k[2,2,2,3]/age.cat[1]

  aging<-p+aging
  aging <<- aging
  
#we need to take into account the larger proportion of the population that is sexually active in the older age group 
#aging.rate is calculated in laod.start, and uses estimates of the size of the not sexually active (NSA) group by age to determine the proportion of the population
#transitioning into the sexually active population as they age from age group 1 to 2. 
#we calculate the proportion of individuals in the NSA group who remain in that category as they age (aging.nsa):
  
  p<-matrix(1,nrow=index, ncol=index) #rate of entry into older age group is split between sexually active and not groups
  p[3:4,1:2]= 1-aging.rate.m[1] #males remaining in NSA group as age
  p[7:8,5:6]=1-aging.rate.m[2]
  p[11:12,9:10]=1-aging.rate.m[3]
  p[15:16,13:14]=(1-aging.rate.m[4])
  p[19:20,13:14]=(1-aging.rate.m[5]) 
  p[19:20,17:18]=1-aging.rate.m[5]
  
  p[23:24,21:22]=1-aging.rate.f[1] #females staying in NSA group as age
  p[27:28,25:26]=1-aging.rate.f[2]
  p[31:32,29:30]=1-aging.rate.f[3]
  
  aging.nsa <<- p  #proportion of NSA individuals remaining in NSA compartment as age group younger to older age group
  
  #rate of entry into youngest age group - we are assuming a constant population size, so rate of births = rate of exit from old age group
  #as older HIV+ MSM exit/reenter the model, they need to be redistributed across HIV- and HIV+ groups to accounted for lower HIV prev in younger age group
  births<-matrix(0,nrow=index, ncol=index) ## rate of births is split between different sized ACs
  births[1,3:4]=p.k[1,1,1,1]/age.cat[2] #brith rate in males
  births[2,3:4]=p.k[1,2,1,1]/age.cat[2]
  births[5,7:8]=p.k[1,1,1,2]/age.cat[2]
  births[6,7:8]=p.k[1,2,1,2]/age.cat[2]
  births[9,11:12]=p.k[1,1,1,3]/age.cat[2]
  births[10,11:12]=p.k[1,2,1,3]/age.cat[2]
  births[13,15:16]=p.k[1,1,1,4]/age.cat[2]
  births[14,15:16]=p.k[1,2,1,4]/age.cat[2]
  births[17,19:20]=(p.k[1,1,1,5]/age.cat[2])*(birth.rate.hiv) #adjusting birth rate of HIV+ MSM to account for lower HIV prev in younger age group
  births[18,19:20]=(p.k[1,2,1,5]/age.cat[2])*(birth.rate.hiv)
  births[13,19:20]=(p.k[1,1,1,4]/age.cat[2])*(1-birth.rate.hiv)
  births[14,19:20]=(p.k[1,2,1,4]/age.cat[2])*(1-birth.rate.hiv)
  
  births[21,23:24]=p.k[2,1,1,1]/age.cat[2] #birth rate in females
  births[22,23:24]=p.k[2,2,1,1]/age.cat[2]
  births[25,27:28]=p.k[2,1,1,2]/age.cat[2]
  births[26,27:28]=p.k[2,2,1,2]/age.cat[2]
  births[29,31:32]=p.k[2,1,1,3]/age.cat[2]
  births[30,31:32]=p.k[2,2,1,3]/age.cat[2]
  births <<- births
  
  # we need to divide these births into the sexually active and NSA groups; p.sa.array is calculated in pop.calc and describes the size of the SA population
  # by age, sex, subpop, and sexual activity group
  # births.sa and birth.nsa are used as multipliers with births to calculate birth rate by sexual activity status
  
  births.sa<-matrix(0,nrow=index, ncol=index) ## rate of births is split between sexually and not sexually active groups
  births.sa[1,3:4]=p.sa.array[1,1,1,1] #proportion in sexually active population, males
  births.sa[2,3:4]=p.sa.array[1,2,1,1]
  births.sa[5,7:8]=p.sa.array[1,1,1,2]
  births.sa[6,7:8]=p.sa.array[1,2,1,2]
  births.sa[9,11:12]=p.sa.array[1,1,1,3]
  births.sa[10,11:12]=p.sa.array[1,2,1,3]
  births.sa[13,15:16]=p.sa.array[1,1,1,4]
  births.sa[14,15:16]=p.sa.array[1,2,1,4]
  births.sa[17,19:20]=p.sa.array[1,1,1,5]
  births.sa[18,19:20]=p.sa.array[1,2,1,5]
  births.sa[13,19:20]=p.sa.array[1,1,1,4] 
  births.sa[14,19:20]=p.sa.array[1,2,1,4] 
  
  births.sa[21,23:24]=p.sa.array[2,1,1,1] #proportion in sexually active population, females
  births.sa[22,23:24]=p.sa.array[2,2,1,1]
  births.sa[25,27:28]=p.sa.array[2,1,1,2]
  births.sa[26,27:28]=p.sa.array[2,2,1,2]
  births.sa[29,31:32]=p.sa.array[2,1,1,3]
  births.sa[30,31:32]=p.sa.array[2,2,1,3]
  births.sa <<- births.sa
  
  
  births.nsa<-matrix(0,nrow=index, ncol=index) 
  births.nsa[1,3:4]=1-p.sa.array[1,1,1,1] # proportion in not sexually active population, males
  births.nsa[2,3:4]=1-p.sa.array[1,2,1,1]
  births.nsa[5,7:8]=1-p.sa.array[1,1,1,2]
  births.nsa[6,7:8]=1-p.sa.array[1,2,1,2]
  births.nsa[9,11:12]=1-p.sa.array[1,1,1,3]
  births.nsa[10,11:12]=1-p.sa.array[1,2,1,3]
  births.nsa[13,15:16]=1-p.sa.array[1,1,1,4]
  births.nsa[14,15:16]=1-p.sa.array[1,2,1,4]
  births.nsa[17,19:20]=1-p.sa.array[1,1,1,5]
  births.nsa[18,19:20]=1-p.sa.array[1,2,1,5]
  births.nsa[13,19:20]=1-p.sa.array[1,1,1,4]
  births.nsa[14,19:20]=1-p.sa.array[1,2,1,4] 
  
  births.nsa[21,23:24]=1-p.sa.array[2,1,1,1]  #proportion in sexually active population, males
  births.nsa[22,23:24]=1-p.sa.array[2,2,1,1]
  births.nsa[25,27:28]=1-p.sa.array[2,1,1,2]
  births.nsa[26,27:28]=1-p.sa.array[2,2,1,2]
  births.nsa[29,31:32]=1-p.sa.array[2,1,1,3]
  births.nsa[30,31:32]=1-p.sa.array[2,2,1,3]
  births.nsa <<- births.nsa
  
} 
