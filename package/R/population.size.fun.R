##################################
### calculate population sizes ###
##################################
#takes total popuation size and calculates the number of individuals in each sex, age, subpopulation, and sexual activity group
#i=subpopulation, j=sex, k=sexual activity group, l=age group

pop.calc <- function(n.total){
  p.s <-c (p.s.1,1-(p.s.1+p.s.3), p.s.3) # proportion of population in subgroups 1, 2, 3
  p.i<- array(0,dim=c(j,k,l,i))  #array of population distribution by subpopulation
  p.i[,,,1] <-p.s[1] #proportion of pop in i=1 (black)
  p.i[,,,3] <-p.s[3] #proportion of pop in i=3 (Hispanic)
  p.i[,,,2] <-p.s[2] #proportion of pop in i=2 (other)
  p.i[1,,1,4] <-p.msm*(1-p.hiv[1]) #proportion of young pop in i=4 (MSM HIV-)
  p.i[1,,2,4] <-p.msm*(1-p.hiv[2]) #proportion of old pop in i=4 (MSM HIV-)
  p.i[1,,1,5]<-p.msm*p.hiv[1] #proportion of young pop in i=5 (MSM HIV+)
  p.i[1,,2,5]<-p.msm*p.hiv[2] #proportion of old pop in i=5 (MSM HIV+)
  p.l<-array(dim=c(j,k,l,i)) #age distribution
  p.l[,,1,]<-age.cat[1]/sum(age.cat) #proportion of population in l=1 (young age cat)
  p.l[,,2,]<-age.cat[2]/sum(age.cat) #proportion of population in l=2 (old age cat)
  p.j<-array(dim=c(j,k,l,i)) #pop dist by sex -- assume equal M and F for het. populations
  p.j[1,,,1:3]<-0.5*(1-p.msm) #male heterosexual population - adjust proportion of heterosexual male population to account for MSM
  p.j[2,,,1:3]<-0.5  #female population
  p.j[1,,,4:5]<-0.5 #MSM population
  p.j[2,,,4:5]<-0 #no females in MSM subpops (4&5)
  p.j.1.1.1 <-rep(p.low,2) #proportion of M and F in AC=1, for i=1, l=1
  p.j.2.2.1 <-rep(p.low,2) #proportion of M and F in AC=1, for i=1 , l=2
  p.j.1.1.2 <-rep(p.low,2) #proportion of M and F in AC=1 for i=2, l=1
  p.j.2.2.2 <-rep(p.low,2) #proportion of M and F in AC=1 for i=2, l=2
  p.j.1.1.3 <-rep(p.low,2) #proportion of M and F in AC=1 for i=3, l=1
  p.j.2.2.3 <-rep(p.low,2) #proportion of M and F in AC=1 for i=3, l=2
  p.j.1.1.4 <-c(p.low.msm, 0) #proportion of M in AC=1 for i=4 (MSM), l=1
  p.j.2.2.4 <-c(p.low.msm,0)  #proportion of M in AC=1 for i=4, l=2
  p.j.1.1.5 <-c(p.low.msm, 0) #proportion of M in AC=1 for i=5 (MSM-HIV+), l=1
  p.j.2.2.5 <-c(p.low.msm,0)  #proportion of M in AC=1 for i=5, l=2
  p.k <- array(c(p.j.1.1.1,  #distribution of population by low and high activity status
                 1-p.j.1.1.1,
                 p.j.2.2.1, 
                 1-p.j.2.2.1, 
                 p.j.1.1.2, 
                 1-p.j.1.1.2, 
                 p.j.2.2.2, 
                 1-p.j.2.2.2,
                 p.j.1.1.3, 
                 1-p.j.1.1.3, 
                 p.j.2.2.3, 
                 1-p.j.2.2.3,
                 p.j.1.1.4,
                 1-p.j.1.1.4,
                 p.j.2.2.4,
                 1-p.j.2.2.4),
               dim=c(j,k,l,i)) #for each i,l, j*k matrix of subpop dist'n
  p.k[2,,,4:5] <- 0 #no females in subpops 4&5 (MSM)
  p.k <<- p.k #used by aging.fun
  
  p.sa.array<-array(dim=c(j,k,l,i)) #proportion of population that is sexually active
  p.sa.array[1,,1,]<-rep(p.sa.m.y, each=2)
  p.sa.array[1,,2,]<-rep(p.sa.m.o, each=2)
  p.sa.array[2,,1,]<-rep(p.sa.f.y, each=2)
  p.sa.array[2,,2,]<-rep(p.sa.f.o, each=2)
  p.sa.array <<- p.sa.array #used by aging.fun

  p.dist.global<-p.i*p.j*p.k*p.l #population distribution across subpopulations
  n.dist<<-n.total*p.dist.global #population size
  n.dist.sa<<-n.dist*p.sa.array #population size for sexually active
  p.dist <<- sweep(p.dist.global,c(3,4),apply(p.dist.global,c(3,4),sum),"/") #population distribution within subpopulation and age group
  p.s.dist <- array(0,dim=c(rep((i),2)))   #relative sizes of different subpopulations - assume that mixing outside of subpopulation is proportionate to number of individuals in each subpopulation, EXCLUDING MSM HERE
  p.s.dist[1,2:3] <-p.s[2:3]/sum(p.s[2:3])
  p.s.dist[2,c(1,3)] <- p.s[c(1,3)]/sum(p.s[c(1,3)])
  p.s.dist[3,1:2] <- p.s[1:2]/sum(p.s[1:2])
  p.s.dist[4,1:3] <- p.s[1:3]/sum(p.s[1:3])
  p.s.dist[5,1:3] <- p.s[1:3]/sum(p.s[1:3])
  p.s.dist <<- p.s.dist
  n.s.dist <<- apply(n.dist,c(1,4),sum)  # population sizes by sex (j) and subpop (i)
  n.s.dist.sa <<- apply(n.dist,c(1,4),sum)  # population sizes by sex (j) and subpop (i) for sexually active population
  n.s.pop <<- c(n.s.dist[1,1:3],n.s.dist[2,1:3]) # pop sizes, minus MSM
  n.s.pop.sa <<- c(n.s.dist.sa[1,1:3],n.s.dist.sa[2,1:3]) # pop sizes, minus MSM, for sexually active population
  n.i<<-c(n.dist[1,,,],n.dist[2,,,]) #(pop size, by sex, subpop, AC, age) M, subpop=1: age=AC=1, age=1 AC=2,age=2 AC=1, age=2 AC=2; M subpop=2: age=AC=1,etc. Then F, subpop1...
  n.sa<<-c(n.dist.sa[1,,,],n.dist.sa[2,,,]) # size of sexually active population
  n.nsa<<- n.i - n.sa # size of not sexually active population
} 