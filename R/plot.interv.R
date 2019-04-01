library(gridExtra)
library(grid)
library(Hmisc)
library(R.utils)

####################################################################
### plot figures comparing the different screening interventions ###
####################################################################
plot.interv.comp <- function(){
pred <- as.data.frame(post.sample$outputs)
#pred <- as.data.frame(t(apply(pred, 1, unlist)))
theta.list<-as.data.frame(post.sample$theta.list)


interv.names = c("2016 levels", "Universal", "Guidelines","Guidelines + Universal", "Prior Infection", "Prior Infection + Universal")
interv.comp <- 1 #intervention to use as comparator when calculating averted cases
interv.comp.name <- "base case"


prev.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=17*n.interv),year=rep(1999:2015, n.interv),interv=rep(unlist(pred$interv),each=17), Male=unlist(pred$prev.m), Female=unlist(pred$prev.f), MSM=unlist(pred$prev.msm)))

prev.end <- subset(prev.out, year==2016)
for(q in 1:max(prev.end$run)) {
  for (w in 1:max(prev.end$interv)) {
    prev.end$pct.avert.m[(q-1)*n.interv +w] <- 100*(prev.end$Male[(q-1)*n.interv + interv.comp] - prev.end$Male[(q-1)*n.interv + w]) / prev.end$Male[(q-1)*n.interv + interv.comp]
    prev.end$pct.avert.f[(q-1)*n.interv +w] <- 100*(prev.end$Female[(q-1)*n.interv + interv.comp] - prev.end$Female[(q-1)*n.interv + w]) / prev.end$Female[(q-1)*n.interv + interv.comp]
    prev.end$pct.avert.msm[(q-1)*n.interv +w] <- 100*(prev.end$MSM[(q-1)*n.interv + interv.comp] - prev.end$MSM[(q-1)*n.interv + w]) / prev.end$MSM[(q-1)*n.interv + interv.comp]
  }
}
prev.end <- melt(prev.end, id.vars=c("year","run","interv" ))
prev.end$Intervention <- factor(prev.end$interv, labels=interv.names)
browser()
# prev.end <- subset(prev.end, (interv!=interv.comp & (variable == "pct.avert.m"| variable=="pct.avert.f" | variable=="pct.avert.msm")), select = -c(year))
# 
# prev.out <- melt(prev.out, id.vars=c("year","run","interv" ))
# prev.out$Intervention <- factor(prev.out$interv, labels=interv.names)
# grouped <- group_by(prev.out, variable, interv, Intervention, year)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE), min=min(value, na.rm=TRUE), max=max(value, na.rm=TRUE))
# 
# prev.out <- subset(prev.out, variable!="MSM")  #plot only M and F
# 
# grouped <- group_by(prev.out, variable, interv, Intervention)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))
# 
# 
# inc.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=17*n.interv) ,year=rep(1999:2015, n.interv), interv=rep(unlist(pred$interv), each=17), Male=unlist(pred$inc.m), Female=unlist(pred$inc.f), MSM=unlist(pred$inc.msm)))
# inc.out <- melt(inc.out, id.vars=c("year","run", "interv"))
# inc.out$Intervention <- factor(inc.out$interv, labels=interv.names)
# 
# grouped <- group_by(inc.out, variable, interv, Intervention, year)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))

#inc.out <- subset(inc.out, variable!="MSM")  #plot only M and F

# diag.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=17*n.interv) ,year=rep(1999:2015, n.interv),interv=rep(unlist(pred$interv), each=17), Male=unlist(pred$diag.m), Female=unlist(pred$diag.f), MSM=unlist(pred$diag.msm)))
# diag.out <- melt(diag.out, id.vars=c("year","run", "interv"))
# diag.out$Intervention <- factor(diag.out$interv, labels=interv.names)
# 
# prev.rr.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015, n.interv), interv=rep(unlist(pred$interv), each=16), "Black M"=unlist(pred$prev.rr.m1), "Black F"=unlist(pred$prev.rr.f1), "Hispanic M"=unlist(pred$prev.rr.m3), "Hispanic F"=unlist(pred$prev.rr.f3),"Other M"=unlist(pred$prev.rr.m2), "Other F"=unlist(pred$prev.rr.f2)))
# prev.rr.out <- melt(prev.rr.out, id.vars=c("year","run", "interv"))
# prev.rr.out$Intervention <- factor(prev.rr.out$interv, labels=interv.names)
# prev.rr.end <- subset(prev.rr.out, year==2015)  #prev rr in 2015
# 
# inc.rr.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015,  n.interv), interv=rep(unlist(pred$interv), each=16), "Black M"=unlist(pred$inc.rr.m1), "Black F"=unlist(pred$inc.rr.f1), "Hispanic M"=unlist(pred$inc.rr.m3), "Hispanic F"=unlist(pred$inc.rr.f3), "Other M"=unlist(pred$inc.rr.m2), "Other F"=unlist(pred$inc.rr.f2)))
# inc.rr.out <- melt(inc.rr.out, id.vars=c("year","run", "interv"))
# inc.rr.out$Intervention <- factor(inc.rr.out$interv, labels=interv.names)
# inc.rr.end <- subset(inc.rr.out, year==2015)  #inc rr in 2015
# grouped <- group_by(inc.rr.end, variable, interv, Intervention)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))
# 
# inc.rr.age.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015,  n.interv), interv=rep(unlist(pred$interv), each=16), "Black My"=unlist(pred$inc.rr.m1y),"Black Mo"=unlist(pred$inc.rr.m1o), "Black Fy"=unlist(pred$inc.rr.f1y), "Black Fo"=unlist(pred$inc.rr.f1o)))
# inc.rr.age.out <- melt(inc.rr.age.out, id.vars=c("year","run", "interv"))
# inc.rr.age.out$Intervention <- factor(inc.rr.age.out$interv, labels=interv.names)
# inc.rr.age.end <- subset(inc.rr.age.out, year==2015)  #inc rr by age in 2015
# grouped <- group_by(inc.rr.age.out, variable, interv, Intervention)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))
# 
# diag.rr.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015, n.interv), interv=rep(unlist(pred$interv), each=16), "Black M"=unlist(pred$diag.rr.m1), "Black F"=unlist(pred$diag.rr.f1), "Hispanic M"=unlist(pred$diag.rr.m3), "Hispanic F"=unlist(pred$diag.rr.f3), "Other M"=unlist(pred$diag.rr.m2), "Other F"=unlist(pred$diag.rr.f2)))
# diag.rr.out <- melt(diag.rr.out, id.vars=c("year","run", "interv"))
# diag.rr.out$Intervention <- factor(diag.rr.out$interv, labels=interv.names)
# diag.rr.end <- subset(diag.rr.out, year==2015)  #diag rr in 2015
# grouped <- group_by(diag.rr.out, variable, interv, Intervention)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))
# 
# diag.rr.age.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015,  n.interv), interv=rep(unlist(pred$interv), each=16), "Black My"=unlist(pred$inc.rr.m1y),"Black Mo"=unlist(pred$inc.rr.m1o), "Black Fy"=unlist(pred$inc.rr.f1y), "Black Fo"=unlist(pred$inc.rr.f1o)))
# diag.rr.age.out <- melt(diag.rr.age.out, id.vars=c("year","run", "interv"))
# diag.rr.age.out$Intervention <- factor(diag.rr.age.out$interv, labels=interv.names)
# diag.rr.age.end <- subset(diag.rr.age.out, year==2015)  #inc rr by age in 2015
# grouped <- group_by(diag.rr.age.out, variable, interv, Intervention)
# inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))

cum.inc.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=1*n.interv) , interv=rep(unlist(pred$interv), each=1), "Male"=unlist(pred$cum.inc.m), "Female"=unlist(pred$cum.inc.f), "MSM"=unlist(pred$cum.inc.msm), "Total"=unlist(pred$cum.inc.tot)))
cum.inc.out$Male <- sapply(pred$prev, function(x) x[["cum.inc.m"]])
cum.inc.out$Female <- sapply(pred$prev, function(x) x[["cum.inc.f"]])
cum.inc.out$MSM <- sapply(pred$prev, function(x) x[["cum.inc.msm"]])
cum.inc.out$Total <- sapply(pred$prev, function(x) x[["cum.inc.tot"]])
# will need this one later
cum.inc.out.ns <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=1*n.interv) , interv=rep(unlist(pred$interv), each=1), "Male"=cum.inc.out$Male, "Female"=cum.inc.out$Female, "MSM"=cum.inc.out$MSM, "Total"=cum.inc.out$Total))

for(q in 1:max(cum.inc.out$run)) {
  for (r in 1:max(cum.inc.out$interv)) {
    cum.inc.out$total.avert[(q-1)*n.interv + r] <-cum.inc.out$Total[(q-1)*n.interv + interv.comp] - cum.inc.out$Total[(q-1)*n.interv + r] 
    cum.inc.out$pct.avert[(q-1)*n.interv + r] <- 100*(cum.inc.out$Total[(q-1)*n.interv + interv.comp] - cum.inc.out$Total[(q-1)*n.interv + r]) / cum.inc.out$Total[(q-1)*n.interv + interv.comp] 
    cum.inc.out$pct.avert.m[(q-1)*n.interv + r] <- 100*(cum.inc.out$Male[(q-1)*n.interv + interv.comp] - cum.inc.out$Male[(q-1)*n.interv + r]) / cum.inc.out$Male[(q-1)*n.interv + interv.comp] 
    cum.inc.out$pct.avert.f[(q-1)*n.interv + r] <- 100*(cum.inc.out$Female[(q-1)*n.interv + interv.comp] - cum.inc.out$Female[(q-1)*n.interv + r]) / cum.inc.out$Female[(q-1)*n.interv + interv.comp] 
    cum.inc.out$pct.avert.msm[(q-1)*n.interv + r] <- 100*(cum.inc.out$MSM[(q-1)*n.interv + interv.comp] - cum.inc.out$MSM[(q-1)*n.interv + r]) / cum.inc.out$MSM[(q-1)*n.interv + interv.comp] 
  }
}
browser()
cum.inc.out <- melt(cum.inc.out, id.vars=c("run", "interv"))
cum.inc.out$Intervention <- factor(cum.inc.out$interv, labels=interv.names)
# cum.inc.avert <- subset(cum.inc.out, ((variable=="pct.avert" | variable=="pct.avert.m" | variable=="pct.avert.f" | variable=="pct.avert.msm" ) & interv!=interv.comp) )

cum.inc.avert <-
  dplyr::filter(
    cum.inc.out,
    interv != interv.comp,
    variable %in% c('pct.avert', 'pct.avert.m', 'pct.avert.f', 'pct.avert.msm')
  )


# cum.inc.avert <- subset(cum.inc.out, interv!=interv.comp, select=c(pct.avert.m, pct.avert.f, pct.avert.msm))
# cum.inc.avert.all <- subset(cum.inc.out, select=c(pct.avert))

cum.inc.avert.all <- filter(cum.inc.avert, variable == 'pct.avert') %>% 
  select(-variable)

#cum.inc.avert.all <- subset(cum.inc.avert.all, select=-variable)
names(cum.inc.avert.all)[names(cum.inc.avert.all)=="value"] <- "cum.inc.avert"
grouped <- group_by(cum.inc.out, interv, variable, Intervention)
#browser()
#inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))

cum.inc.avert$value <- as.numeric(cum.inc.avert$value)
scaleFUN <- function(x) sprintf("%.2f", x)
plot.cum.inc <- ggplot(data=cum.inc.avert, aes(x=variable, y=value, group=interaction(Intervention, variable)))+
 stat_boxplot(aes(fill=Intervention), outlier.shape=NA) +
    labs(x="", y="Incident cases averted (%)\n(relative to screening at 2016 levels)")+
    #coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(breaks = seq(-300, 100, len = 50)) +
    theme_bw()+
  theme(legend.position="right",  axis.text.x=element_text(angle=0,hjust=0.5,size=14),axis.text.y=element_text(size=14),  
       plot.title=element_text(size=14, hjust=0),
       axis.title = element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
       scale_x_discrete(labels=c("pct.avert" = "Total", "pct.avert.m" = "Males", "pct.avert.f" = "Females",
                            "pct.avert.msm" = "MSM only")) +
      scale_y_continuous(labels=scaleFUN)
browser()
pct.avert <- rbind(prev.end, cum.inc.avert)
#grouped <- group_by(pct.avert, variable, interv, Intervention)
#inc.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))

p.msm.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2011:2016, n.interv), interv=rep(unlist(pred$interv), each=16), "Prevalent cases"=unlist(pred$p.prev.msm), "Incident cases"=unlist(pred$p.inc.msm), "Reported cases"=unlist(pred$p.diag.msm)))
p.msm.end <- subset(p.msm.out, year==2016)  #diag rr in 2015
for(q in 1:max(p.msm.end$run)) {
  for (r in 1:max(p.msm.end$interv)) {
    p.msm.end$inc.pct.avert[(q-1)*n.interv + r] <- 100*(p.msm.end$`Incident cases`[(q-1)*n.interv + r] - p.msm.end$`Incident cases`[(q-1)*n.interv + interv.comp] ) / p.msm.end$`Incident cases`[(q-1)*n.interv + interv.comp]
  }
}

p.msm.out <- melt(p.msm.out, id.vars=c("year","run", "interv"))
p.msm.out$Intervention <- factor(p.msm.out$interv, labels=interv.names)
p.msm.end <- melt(p.msm.end, id.vars=c("year","run", "interv"))
p.msm.end$Intervention <- factor(p.msm.end$interv, labels=interv.names)
grouped <- group_by(p.msm.end, variable, interv, Intervention)
msm.sum <- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))


p.msm.end.inc <- subset(p.msm.end[,2:6], variable=="Incident cases")
p.msm.end.inc <- subset(p.msm.end.inc, select=-variable)
names(p.msm.end.inc)[names(p.msm.end.inc)=="value"] <- "p.msm"
msm.inc<-merge(p.msm.end.inc, cum.inc.avert.all)

p.msm.end.inc.pct <- subset(p.msm.end[,2:6], variable=="inc.pct.avert")
p.msm.end.inc.pct <- subset(p.msm.end.inc.pct, select=-variable)
names(p.msm.end.inc.pct)[names(p.msm.end.inc.pct)=="value"] <- "p.msm.pct"
msm.inc.pct<-merge(p.msm.end.inc.pct, cum.inc.avert.all)

test.vol <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=n.interv) , interv=rep(unlist(pred$interv), each=1), "tests"=unlist(pred$test.tot)))



#compare to no screening for number needed to screen
#cum.inc.out.ns <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=1*n.interv) , interv=rep(unlist(pred$interv), each=1), "Male"=unlist(pred$cum.inc.m), "Female"=unlist(pred$cum.inc.f), "MSM"=unlist(pred$cum.inc.msm), "Total"=unlist(pred$cum.inc.tot)))
# test <- list()
# for(q in 1:max(cum.inc.out.ns$run)) {
#   for (r in 1:max(cum.inc.out.ns$interv)) {
#     test$total.avert[(q-1)*n.interv + r] <-cum.inc.out.ns$Total[(q-1)*n.interv + 1] - cum.inc.out.ns$Total[(q-1)*n.interv + r]
#     test$pct.avert[(q-1)*n.interv + r] <- 100*(cum.inc.out.ns$Total[(q-1)*n.interv + 1] - cum.inc.out.ns$Total[(q-1)*n.interv + r]) / cum.inc.out$Total[(q-1)*n.interv + 1]
#   }
# }
# 
# cum.inc.out.ns <- melt(cum.inc.out.ns, id.vars=c("run", "interv"))
# cum.inc.out.ns$Intervention <- factor(cum.inc.out.ns$interv, labels=interv.names)
# 
# 
# cum.cases.avert.ns <- subset(cum.inc.out.ns, variable=="total.avert")
# nnt <- merge(test.vol, cum.cases.avert.ns, by=c("interv","run"))
# nnt$nnt <- nnt$tests/nnt$value
# nnt <- subset(nnt, interv!=1)
# 
# grouped <- group_by(nnt, variable, interv, Intervention)
# nnt.sum<- summarise(grouped, mean=mean(value, na.rm=TRUE), p2.5=quantile(value,0.025, na.rm=TRUE), p97.5=quantile(value,0.975, na.rm=TRUE))
# 
# 
# inc.age.subpop <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=16*n.interv) ,year=rep(2000:2015, n.interv), interv=rep(unlist(pred$interv), each=16), "Black M 20-44y"=unlist(pred$inc.m1y), "Black M 45-64y"=unlist(pred$inc.m1o),
#                                       "Other M 20-44y"=unlist(pred$inc.m2y), "Other M 45-64y"=unlist(pred$inc.m2o),
#                                       "Hispanic M 20-44y"=unlist(pred$inc.m3y),"Hispanic M 45-64y"=unlist(pred$inc.m3o),
#                                       "Black F 20-44y"=unlist(pred$inc.f1y), "Black F 45-64y"=unlist(pred$inc.f1o),
#                                       "Other F 20-44y"=unlist(pred$inc.f2y), "Other F 45-64y"=unlist(pred$inc.f2o),
#                                       "Hispanic F 20-44y"=unlist(pred$inc.f3y), "Hispanic F 45-64y"=unlist(pred$inc.f3o)))
# 
# 
# inc.age.out <- as.data.frame(cbind(run=rep(1:max(unlist(pred$iter)),each=17*n.interv) ,year=rep(1999:2015, n.interv), interv=rep(unlist(pred$interv), each=17), "Male 20-44y"=unlist(pred$inc.y.m), "Female 20-44y"=unlist(pred$inc.y.f), "Male 45-64y"=unlist(pred$inc.o.m), "Female 45-64y"=unlist(pred$inc.o.f)))
# inc.age.out <- melt(inc.age.out, id.vars=c("year","run", "interv"))
# inc.age.out$Intervention <- factor(inc.age.out$interv, labels=interv.names)
# inc.age.out<-cbind(inc.age.out, colsplit(inc.age.out$variable, " ", names = c('Sex', 'Age')))
# grouped <- group_by(inc.age.out, variable, interv, Intervention, year)
# #inc.sum<- summarise(grouped, mean=mean(value), p2.5=quantile(value,0.025), p97.5=quantile(value,0.975))
# 
# 
# ### print annual incidence ###
# inc.subpop<-melt(inc.age.subpop, id.vars=c("year", "run", "interv"))
# inc.subpop$value <- signif(inc.subpop$value,3)
# inc.table<-aggregate(inc.subpop$value, FUN=mean, na.rm=TRUE,
#                      by= list(year=inc.subpop$year, interv=inc.subpop$interv, pop=inc.subpop$variable))
# inc.table <- subset(inc.table, interv==3)
# inc.tab <- dcast(inc.table, year~pop, value.var="x")
# write.csv(inc.tab, file="Model_outputs/inc_table.csv")
#####################################################

### Figure 3A ###
# inc.out$Sex <- relevel(inc.out$variable, "Female")
# plot.inc.best <- ggplot(data=subset(inc.out, (interv==3 & year!=1999)), aes(x=year, y=value*100, group=Sex)) +
#   stat_summary(aes(colour=Sex,fill=Sex), geom="smooth",
#                fun.y="mean",
#                fun.ymin = function(z) {quantile(z,0.025)},
#                fun.ymax = function(z) {quantile(z,0.975)}) +
#   coord_cartesian(ylim=c(0,4)) +
#   labs(x="Year", y="Syphilis incidence (%)\n\n") +
#   theme_classic() +
#   theme(legend.title=element_blank(), axis.text.x=element_text(angle=0,hjust=0.5,size=10),axis.text.y=element_text(size=10),
#         plot.title=element_text(size=16, hjust=0),
#         legend.justification = c(1, 1), legend.position = c(1, 1),
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
#   ggtitle("A")
# 
# ### Figure 3B ###
# plot.inc.age.1 <- ggplot(data=subset(inc.age.out, (interv==3 & year!=1999)), aes(x=year, y=value*100, group=variable)) +
#   stat_summary(aes(colour=Age,fill=Age), geom="smooth",
#                fun.y="mean",
#                fun.ymin = function(z) {quantile(z,0.025)},
#                fun.ymax = function(z) {quantile(z,0.975)}) +
#   facet_grid(.~Sex) +
#   coord_cartesian(ylim=c(0,4)) +
#   labs(x="Year", y="Syphilis incidence (%)\n") +
#   theme_classic() +
#   theme(legend.position="right", legend.title=element_blank(), axis.text.x=element_text(angle=0,hjust=0.5,size=10),axis.text.y=element_text(size=10),
#         plot.title=element_text(size=16, hjust=0),
#         panel.spacing = unit(1, "lines"),
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
#   scale_fill_manual(values=cols[6:7])+
#   scale_color_manual(values=cols[6:7])+
#   ggtitle("B")
# 
# ### Figure 3C ###
# plot.p.msm.best <- ggplot(data=subset(p.msm.out, (interv==3 & variable=="Incident cases")), aes(x=year, y=value)) +
#   stat_summary( aes(colour=variable, fill=variable), geom="smooth",
#                 fun.y="mean",
#                 fun.ymin = function(z) {quantile(z,0.025)},
#                 fun.ymax = function(z) {quantile(z,0.975)}) +
#   coord_cartesian(ylim=c(0,1)) +
#   labs(x="Year", y="Proportion of incident male cases in MSM\n") +
#   theme_classic()+
#   theme(legend.position="none",  axis.text.x=element_text(angle=0,hjust=0.5,size=10),axis.text.y=element_text(size=10),
#         plot.title=element_text(size=16, hjust=0),
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
#   ggtitle("C")
# 
# ### Figure 3D ###
# plot.inc.rr.best <- ggplot(data=subset(inc.rr.out, interv==3), aes(x=year, y=value, group=variable)) +
#   stat_summary(aes(colour=variable,fill=variable), geom="smooth",
#                fun.y="mean",
#                fun.ymin = function(z) {quantile(z,0.025)},
#                fun.ymax = function(z) {quantile(z,0.975)}) +
#   expand_limits(y=0) +
#   labs(x="Year", y="Incidence rate ratio\n") +
#   theme_classic() +
#   theme(legend.title=element_blank(), axis.text.x=element_text(angle=0,hjust=0.5,size=10),axis.text.y=element_text(size=10),
#         plot.title=element_text(size=16, hjust=0),
#         legend.position = "right",
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
#   ggtitle("D")


### Figure 4 ###
msm.sub <- subset(msm.inc.pct, (interv!=1 & interv!=2))
msm.sub.mean <- msm.sub %>%
  group_by(Intervention) %>%
  summarise(p.msm.pct = mean(p.msm.pct, na.rm=TRUE),
           cum.inc.avert = mean(cum.inc.avert, na.rm=TRUE))

plot.p.msm.inc.pct <- ggplot(data=subset(msm.inc.pct, (interv!=1 & interv!=2)), aes(x=cum.inc.avert, y=p.msm.pct, color=Intervention)) +
  geom_rect(aes(xmin = 0, xmax = 100, ymin = -Inf, ymax = 0), alpha=0.1, fill="grey95", color = NA) +
  geom_point(size=1, shape=16, alpha=0.3) +
  geom_point(data=msm.sub.mean, shape=3, size=2, stroke=3) +
  expand_limits(y=0) +
  labs(x=bquote(atop("Incident cases averted (%)", "(relative to"~ .(interv.comp.name)*")")),
       y=bquote(atop("Change in male cases occurring in MSM (%)", "(relative to"~ .(interv.comp.name)*")")) )+
  theme_classic()+
  theme(legend.position="right",  axis.text.x=element_text(angle=0,hjust=1,size=12),axis.text.y=element_text(size=12),
        plot.title=element_text(size=16, hjust=0),
        axis.title = element_text(size=16),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  guides(color=guide_legend(override.aes = list(shape=16, size=1, alpha=0.7, stroke=1)))


### Figure 5A
plot.inc.rr.end <- ggplot(data=subset(inc.rr.end, (interv!=1 & interv!=2)), aes(x=variable, y=value, group=interaction(Intervention, variable)))+
  stat_boxplot(aes(fill=Intervention), outlier.shape=NA) +
  labs(x="\nPopulation", y="Incidence rate ratio in 2016\n")+
  geom_hline(yintercept=1,color="grey70",lty=3)+
  coord_cartesian(ylim=c(0,5))+
  theme_classic()+
  #coord_flip()+
  theme(legend.position="right",  axis.text.x=element_text(angle=0,hjust=0.5,size=14),axis.text.y=element_text(size=14),
        plot.title=element_text(size=16, hjust=0),
        axis.title = element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  ggtitle("A")

legend.rr <- get_legend(plot.inc.rr.end)
plot.inc.rr.end <- plot.inc.rr.end + theme(legend.position ="none")

### Figure 5B
plot.diag.rr.end <- ggplot(data=subset(diag.rr.end, (interv!=1 & interv!=2)), aes(x=variable, y=value, group=interaction(Intervention, variable)))+
  stat_boxplot(aes(fill=Intervention), outlier.shape=NA) +
  labs(x="\nPopulation", y="Reported case rate ratio in 2016\n") +
  geom_hline(yintercept=1,color="grey70",lty=3)+
  coord_cartesian(ylim=c(0,5))+
  theme_classic()+
  theme(legend.position="none",  axis.text.x=element_text(angle=0,hjust=0.5,size=14),axis.text.y=element_text(size=14),
        plot.title=element_text(size=16, hjust=0),
        axis.title = element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  ggtitle("B")

### Figure 6 ###

nnt <- subset(nnt, (interv!=1 & interv!=2))
plot.nns <- ggplot(data=nnt, aes(x=Intervention, y=nnt))+
  stat_boxplot(aes(fill=Intervention), outlier.shape=NA) +
  labs(x="", y="Number needed to screen to avert one syphilis infection\n") +
  expand_limits(y=1) +
  theme_classic()+
  theme(legend.position="none",  axis.text.x=element_text(angle=0,hjust=0.5,size=14),axis.text.y=element_text(size=14),
        axis.title = element_text(size=14),
        plot.title=element_text(size=16, hjust=0),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  coord_flip()

# ####################
# ### save figures ###
# ####################
#browser()
png('Model_outputs/interv_eg.png', units="in", width=8, height=6, res=600)
plot.cum.inc
dev.off()

png('Model_outputs/F4.png', units="in", width=10, height=10, res=300)
grid.arrange(plot.p.msm.inc.pct,nrow=1)
dev.off()

png('Model_outputs/F5.png', units="in", width=10, height=10, res=300)
grid.arrange(plot.inc.rr.end, legend.rr, plot.diag.rr.end, widths=c(3,1), ncol=2)
dev.off()

png('Model_outputs/F6.png', units="in", width=10, height=10, res=300)
grid.arrange(plot.nns,nrow=1)
dev.off()

}