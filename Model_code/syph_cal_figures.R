
pred <- as.data.frame(post.sample$outputs)
diag.len <- nrow(diag.age.sex.dat)
msm.len <- nrow(msm.dat) #years of data for proportion of male cases reported as msm
hiv.len <- nrow(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"])))) #years of data for hiv coinfection
out.subpop <- melt(as.matrix(subset(pred, select=subpop.assort1:subpop.assort8))) #model subpopulation assortativity 
out.age <- melt(as.matrix(subset(pred, select=age.assort1:age.assort4))) #model age assortativity
out.diag.y.m <- melt(as.matrix(subset(pred, select=prev.early.inf.rate.m1:(prev.early.inf.rate.m1+diag.len-1))))  #reported early syph cases by age-sex
out.diag.o.m <- melt(as.matrix(subset(pred, select=(prev.early.inf.rate.m1+diag.len):(prev.early.inf.rate.m1+diag.len*2-1))))  
out.diag.y.f <- melt(as.matrix(subset(pred, select=prev.early.inf.rate.f.y1:(prev.early.inf.rate.f.y1+diag.len-1)))) 
out.diag.o.f <- melt(as.matrix(subset(pred, select=prev.early.inf.rate.f.o1:(prev.early.inf.rate.f.o1+diag.len-1)))) 
out.diaglate.y.m <- melt(as.matrix(subset(pred, select=prev.diag.late1:(prev.diag.late1+diag.len-1))))  #reported late syph cases by age-sex
out.diaglate.o.m <- melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len):(prev.diag.late1+diag.len*2-1))))  
out.diaglate.y.f <- melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len*2):(prev.diag.late1+diag.len*3-1))))
out.diaglate.o.f <- melt(as.matrix(subset(pred, select=(prev.diag.late1+diag.len*3):(prev.diag.late1+diag.len*4-1))))  
syph.ratio.y.m <- out.diag.y.m$value/(out.diaglate.y.m$value+ out.diag.y.m$value) #proportion of early cases to all cases by age-sex
syph.ratio.o.m <- out.diag.o.m$value/(out.diaglate.o.m$value+ out.diag.o.m$value)
syph.ratio.y.f <- out.diag.y.f$value/(out.diaglate.y.f$value+ out.diag.y.f$value)
syph.ratio.o.f <- out.diag.o.f$value/(out.diaglate.o.f$value+ out.diag.o.f$value)

out.rr.diag <- melt(as.matrix(subset(pred, select=prev.diag.rr1:prev.diag.rr8))) ## reported case relative risk (pooled estimates for last 5 years)
out.p.msm.y <- melt(as.matrix(subset(pred, select=prev.p.diag.msm1:(prev.p.diag.msm1+msm.len-1)))) #proportion of male cases in young MSM
out.p.msm.o <- melt(as.matrix(subset(pred, select=(prev.p.diag.msm1+msm.len):(prev.p.diag.msm1+msm.len*2-1)))) #proportion of male cases in old MSM
out.p.hiv.y <- melt(as.matrix(subset(pred, select=prev.p.diag.hiv1:(prev.p.diag.hiv1+hiv.len-1)))) #proportion of young MSM cases with HIV
out.p.hiv.o <- melt(as.matrix(subset(pred, select=(prev.p.diag.hiv1+hiv.len):(prev.p.diag.hiv1+hiv.len*2-1)))) #proportion of old MSM cases with HIV
out.p.sec.y.m <- melt(as.matrix(subset(pred, select=prev.p.diag.sec1:(prev.p.diag.sec1+diag.len-1)))) #proportion of early cases that are secondary
out.p.sec.o.m <- melt(as.matrix(subset(pred, select=(prev.p.diag.sec1+diag.len):(prev.p.diag.sec1+diag.len*2-1)))) #proportion of early cases that are secondary
out.p.sec.y.f <- melt(as.matrix(subset(pred, select=(prev.p.diag.sec1+diag.len*2):(prev.p.diag.sec1+diag.len*3-1)))) #proportion of early cases that are secondary
out.p.el.y.m <- melt(as.matrix(subset(pred, select=prev.p.diag.el1:(prev.p.diag.el1+diag.len-1)))) #proportion of early cases that are early latent
out.p.el.o.m <- melt(as.matrix(subset(pred, select=(prev.p.diag.el1+diag.len):(prev.p.diag.el1+diag.len*2-1)))) #proportion of early cases that are early latent
out.p.el.y.f <- melt(as.matrix(subset(pred, select=(prev.p.diag.el1+diag.len*2):(prev.p.diag.el1+diag.len*3-1)))) #proportion of early cases that are early latent
out.p.early <- melt(as.matrix(subset(pred, select=(prev.p.diag.early)))) #proportion of early cases that are early latent

out.inc.y.m.1 <-melt(100*as.matrix(subset(pred, select=prev.n.inc1:(prev.n.inc1+cal.period-1)))/sum(n.i[y.m.1])) #model incidence by age, sex, and subpopulation
out.inc.o.m.1 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period):(prev.n.inc1+cal.period*2-1)))/sum(n.i[o.m.1]))
out.inc.y.m.2 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*2):(prev.n.inc1+cal.period*3-1)))/sum(n.i[y.m.2]))
out.inc.o.m.2 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*3):(prev.n.inc1+cal.period*4-1)))/sum(n.i[o.m.2]))
out.inc.y.m.3 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*4):(prev.n.inc1+cal.period*5-1)))/sum(n.i[y.m.3]))
out.inc.o.m.3 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*5):(prev.n.inc1+cal.period*6-1)))/sum(n.i[o.m.3]))
out.inc.y.msm <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*6):(prev.n.inc1+cal.period*7-1)))/sum(n.i[y.m.4]))
out.inc.o.msm <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*7):(prev.n.inc1+cal.period*8-1)))/sum(n.i[o.m.4]))
out.inc.y.msmhiv <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*8):(prev.n.inc1+cal.period*9-1)))/sum(n.i[y.m.5]))
out.inc.o.msmhiv <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*9):(prev.n.inc1+cal.period*10-1)))/sum(n.i[o.m.5]))
out.inc.y.f.1 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*10):(prev.n.inc1+cal.period*11-1)))/sum(n.i[y.f.1]))
out.inc.o.f.1 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*11):(prev.n.inc1+cal.period*12-1)))/sum(n.i[o.f.1]))
out.inc.y.f.2 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*12):(prev.n.inc1+cal.period*13-1)))/sum(n.i[y.f.2]))
out.inc.o.f.2 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*13):(prev.n.inc1+cal.period*14-1)))/sum(n.i[o.f.2]))
out.inc.y.f.3 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*14):(prev.n.inc1+cal.period*15-1)))/sum(n.i[y.f.3]))
out.inc.o.f.3 <-melt(100*as.matrix(subset(pred, select=(prev.n.inc1+cal.period*15):(prev.n.inc1+cal.period*16-1)))/sum(n.i[o.f.3]))

#prep input data for plotting
diag.data <- as.data.frame(cbind(year=seq((end.year-(nrow(diag.age.sex.dat))+1),end.year,1), rate=as.numeric(c(diag.age.sex.rate.m, diag.age.sex.rate.f.y, diag.age.sex.rate.f.o)), cat=rep(c("y.m", "o.m", "y.f", "o.f"),each=nrow(diag.age.sex.dat))))
age.data <- as.data.frame(age.dist.dat[1:4,])

#get max values of reported cases to set y-axis
max.y.m <- max(out.diag.y.m$value, na.rm=TRUE)
max.o.m <- max(out.diag.o.m$value, na.rm=TRUE)
max.y.f <-max(out.diag.y.f$value, na.rm=TRUE)
max.o.f <-max(out.diag.o.f$value, na.rm=TRUE)
max.diag <- max(max.y.m, max.o.m, max.y.f, max.o.f,diag.data$diag.age.sex.rate) +5

### plot reported cases by age and sex ###
plot.diag.y.m <- ggplot(data=out.diag.y.m)+
  geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
  stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
  theme_classic() + 
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_m), color="red", shape=15, size=2) +
  labs(title="All M 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
  coord_cartesian(ylim=c(0,max.diag))+
  #expand_limits(y=0)+
  scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))

plot.diag.o.m <- ggplot(data=out.diag.o.m)+
  geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
  stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
  theme_classic() + 
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_m), color="red", shape=15, size=2) +
  labs(title="All M 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
  coord_cartesian(ylim=c(0,max.diag))+
  #expand_limits(y=0)+
  scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))

plot.diag.y.f <- ggplot(data=out.diag.y.f)+
  geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
  stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
  theme_classic() + 
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_y_f), color="red", shape=15, size=2) +
  labs(title="All F 20-44 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
  coord_cartesian(ylim=c(0,max.diag))+
  #expand_limits(y=0)+
  scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))

plot.diag.o.f <- ggplot(data=out.diag.o.f)+
  geom_line(aes(x=Var2, y=value, group=Var1), color="grey") +
  stat_summary(aes(x=Var2, y=value, group=1), geom="line", fun.y="median", color="black", size=1) +
  theme_classic() + 
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  geom_point(data=as.data.frame(diag.age.sex.dat),aes(x=1:diag.len, y=diag_all_o_f), color="red", shape=15, size=2) +
  labs(title="All F 45-64 y", x="Year", y="Reported early syphilis cases\nper 100,000") +
  coord_cartesian(ylim=c(0,max.diag))+
  ##expand_limits(y=0)+
  scale_x_discrete(labels=seq((end.year-diag.len+1),end.year,1))

### plot MSM and HIV data ###
plot.p.msm.y <- ggplot(data=out.p.msm.y)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_point(data=as.data.frame(p.msm.cases[1:msm.len]),aes(x=1:msm.len, y=p.msm.cases[1:msm.len]), color="red", shape=15, size=3) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="MSM 20-44 y", x="Year", y="Proportion of early syphilis\ncases in MSM") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-msm.len+1,end.year,1))

plot.p.msm.o <- ggplot(data=out.p.msm.o)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_point(data=as.data.frame(p.msm.cases[(msm.len+1):(2*msm.len)]),aes(x=1:msm.len, y=p.msm.cases[(msm.len+1):(2*msm.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="MSM 45-64 y", x="Year", y="Proportion of early syphilis\ncases in MSM") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-msm.len+1,end.year,1))

plot.p.hiv.y <- ggplot(data=out.p.hiv.y)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_point(data=as.data.frame(p.hiv.cases[1:hiv.len]),aes(x=1:hiv.len, y=p.hiv.cases[1:hiv.len]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="MSM 20-44 y", x="Year", y="Proportion of MSM cases\nwith HIV coinfection") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-hiv.len+1,end.year,1))

plot.p.hiv.o <- ggplot(data=out.p.hiv.o)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_point(data=as.data.frame(p.hiv.cases[(hiv.len+1):(2*hiv.len)]),aes(x=1:hiv.len, y=p.hiv.cases[(hiv.len+1):(2*hiv.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="MSM 45-64 y", x="Year", y="Proportion of MSM cases\nwith HIV coinfection") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-hiv.len+1,end.year,1))

plot.p.sec.y.m <- ggplot(data=out.p.sec.y.m)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.sec[1:diag.len]),aes(x=1:diag.len, y=p.sec[1:diag.len]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="M 20-44 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.sec.o.m <- ggplot(data=out.p.sec.o.m)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.sec[(diag.len+1):(2*diag.len)]),aes(x=1:diag.len, y=p.sec[(diag.len+1):(2*diag.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=12), axis.text.x=element_text(angle=90,hjust=1, size=12), title=element_text(size=12)) + 
  labs(title="M 45-64 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.sec.y.f <- ggplot(data=out.p.sec.y.f)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.sec[(2*diag.len+1):(3*diag.len)]),aes(x=1:diag.len, y=p.sec[(2*diag.len+1):(3*diag.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=12), axis.text.x=element_text(angle=90,hjust=1, size=12), title=element_text(size=12)) + 
  labs(title="F 20-44 y", x="Year", y="Proportion of reported early cases\nthat are secondary") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.el.y.m <- ggplot(data=out.p.el.y.m)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.el[1:diag.len]),aes(x=1:diag.len, y=p.el[1:diag.len]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(title="M 20-44 y", x="Year", y="Proportion of reported early cases\nthat are early latent") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.el.o.m <- ggplot(data=out.p.el.o.m)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.el[(diag.len+1):(2*diag.len)]),aes(x=1:diag.len, y=p.el[(diag.len+1):(2*diag.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=12), axis.text.x=element_text(angle=90,hjust=1, size=12), title=element_text(size=12)) + 
  labs(title="M 45-64 y", x="Year", y="Proportion of diagnosed early cases\nthat are early latent") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.el.y.f <- ggplot(data=out.p.el.y.f)+
  geom_boxplot(aes(x=Var2, y=value), fill="purple1") +
  geom_point(data=as.data.frame(p.el[(2*diag.len+1):(3*diag.len)]),aes(x=1:diag.len, y=p.el[(2*diag.len+1):(3*diag.len)]), color="red", shape=15, size=2) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=12), axis.text.x=element_text(angle=90,hjust=1, size=12), title=element_text(size=12)) + 
  labs(title="F 20-44 y", x="Year", y="Proportion of diagnosed early cases\nthat are early latent") +
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(labels=seq(end.year-diag.len+1,end.year,1))

plot.p.early <- ggplot(data=out.p.early)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_pointrange(data=p.early.dat,aes(x=1, y=mean, ymin=min, ymax=max), color="red", shape=15, size=0.5) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs(x="", y="Proportion of all reported syphilis cases \nthat are primary, secondary, or early latent") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_x_discrete(labels="")

plot.rr.diag <- ggplot(data=out.rr.diag)+
  geom_boxplot(aes(x=Var2, y=value, fill=Var2)) +
  geom_pointrange(data=diag.rr, aes(x=1:8, y=mean, ymin=min, ymax=max), color="red", shape=15, size=0.75, lty=3) +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=16), axis.text.x=element_text(angle=90,hjust=1, size=16), title=element_text(size=16)) + 
  labs( x="Population", y="Reported case relative risk") +
  scale_x_discrete(labels=c("M black y","M black o", "M Hispanic y", "M Hispanic o", "F black y", "F black o", "F Hispanic y", "F Hispanic o"))


### save plots to pdf in Model_outputs folder ###
png(file="Model_outputs/syph_cal_diag_MA.png", units="in", width=12, height=8, res=600)
grid.arrange(plot.diag.y.m, plot.diag.o.m,
             plot.diag.y.f, plot.diag.o.f,
             plot.p.msm.y, plot.p.msm.o,
             plot.p.hiv.y, plot.p.hiv.o,
             ncol=4)
dev.off()

png(file="Model_outputs/syph_cal_trt_MA.png", units="in", width=8, height=8, res=600)
grid.arrange(plot.p.sec.y.m,plot.p.el.y.m, 
             plot.p.early, plot.rr.diag, 
             ncol=2)
dev.off()


