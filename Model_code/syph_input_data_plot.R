##########################################
###   view syphilis model input data   ###
###   for Louisiana and Massachusetts  ###
##########################################
# update working directory before running
# make sure functions are saved in 'Model_code' folder in directory
# make sure data inputs are saved in 'R_inputs' folder in directory

setwd("~/Dropbox/PPML syphilis/")
source("Model_code/load.start.conditions.R")


############################################
### load required packages and functions ###
############################################
load.start()

library(gridExtra)
library(grid)
library(Hmisc)
library(R.utils)
bezier <- bezier::bezier 

get_legend<-function(myggplot){   ### function to get legend from plot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blankPlot <- ggplot()+geom_blank(aes(1,1)) +  #blank plot to be used as a place holder
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )


### ggplot colors ###
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 10
cols = gg_color_hue(n)

###############################################
### plot data for LA and MA and save to pdf ###
###############################################

states <- c("LA","MA")  #LA= Louisiana; MA= Massachusetts

for (i in 1:length(states)){
state <- states[i] 
statename <- if(state=="LA") "Louisiana" else "Massachusetts"
load.start()
bezier <- bezier::bezier 

diag.age.sex <- tibble::rownames_to_column(as.data.frame(diag.age.sex.dat), "year")
diag.dat <- melt(diag.age.sex[,1:5],  id.vars = "year")
diag.dat$cat <- ifelse(diag.dat$variable=="diag_all_y_m", "M 20-44 y", ifelse(diag.dat$variable=="diag_all_o_m", "M 45-64 y", ifelse(diag.dat$variable=="diag_all_y_f", "F 20-44 y", "F 45-64 y")))
pmsm.dat <- tibble::rownames_to_column(as.data.frame(msm.dat[,c("pMSM_y","pMSM_o")]), "year")
pmsm.dat <- melt(pmsm.dat,  id.vars = "year")
pmsm.dat$cat <- ifelse(pmsm.dat$variable=="pMSM_y", "M 20-44 y", "M 45-64 y" ) 
hiv.dat <- tibble::rownames_to_column(as.data.frame(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"])))), "year") 
hiv.dat <- melt(hiv.dat,  id.vars = "year")
hiv.dat$cat <- ifelse(hiv.dat$variable=="pHIV_y", "M 20-44 y", "M 45-64 y" ) 
sec.dat <- tibble::rownames_to_column(as.data.frame(stage.dat[,c("pSec_y_m","pSec_o_m","pSec_y_f", "pSec_o_f")]), "year" )#proportion cases diagnosed with secondary syphilis
sec.dat <- melt(sec.dat,  id.vars = "year")
sec.dat$cat <- ifelse(sec.dat$variable=="pSec_y_m", "M 20-44 y", ifelse(sec.dat$variable=="pSec_o_m", "M 45-64 y", ifelse(sec.dat$variable=="pSec_y_f", "F 20-44 y", "F 45-64 y")))
el.dat  <- tibble::rownames_to_column(as.data.frame(stage.dat[,c("pEL_y_m","pEL_o_m","pEL_y_f", "pEL_o_f")]), "year") #proportion cases diagnosed with early latent syphilis
el.dat <- melt(el.dat,  id.vars = "year")
el.dat$cat <- ifelse(el.dat$variable=="pEL_y_m", "M 20-44 y", ifelse(el.dat$variable=="pEL_o_m", "M 45-64 y", ifelse(el.dat$variable=="pEL_y_f", "F 20-44 y", "F 45-64 y")))
subpop.dat <- melt(diag.subpop, id.vars="year")
subpop.dat$cat <- ifelse(subpop.dat$variable=="y_m_b", "Black M 20-44 y", 
                         ifelse(subpop.dat$variable=="y_m_o", "Other M 20-44 y", 
                                ifelse(subpop.dat$variable=="y_m_h", "Hispanic M 20-44 y", 
                                       ifelse(subpop.dat$variable=="o_m_b", "Black M 45-64 y", 
                                              ifelse(subpop.dat$variable=="o_m_o", "Other M 45-64 y", 
                                                     ifelse(subpop.dat$variable=="o_m_h", "Hispanic M 45-64 y", 
                                                            ifelse(subpop.dat$variable=="y_f_b", "Black F 20-44 y", 
                                                                   ifelse(subpop.dat$variable=="y_f_o", "Other F 20-44 y", 
                                                                          ifelse(subpop.dat$variable=="y_f_h", "Hispanic F 20-44 y", 
                                                                                 ifelse(subpop.dat$variable=="o_f_b", "Black F 45-64 y",
                                                                                        ifelse(subpop.dat$variable=="o_f_o", "Other F 45-64 y", "Hispanic F 45-64 y")))))))))))
subpop.dat$agecat <- ifelse(subpop.dat$variable=="y_m_b"|subpop.dat$variable=="y_m_o"|subpop.dat$variable=="y_m_h"|
                            subpop.dat$variable=="y_f_b"|subpop.dat$variable=="y_f_o"|subpop.dat$variable=="y_f_h", "20-44 y", "45-64 y" )
                                                            
subpop.dat$sex <- ifelse(subpop.dat$variable=="y_m_b"|subpop.dat$variable=="y_m_o"|subpop.dat$variable=="y_m_h"|
                         subpop.dat$variable=="o_m_b"|subpop.dat$variable=="o_m_o"|subpop.dat$variable=="o_m_h", "M", "F" )

subpop.dat$subpop <- ifelse(subpop.dat$variable=="y_m_b"|subpop.dat$variable=="o_m_b"|subpop.dat$variable=="y_f_b"|subpop.dat$variable=="o_f_b", "Black", 
                            ifelse(subpop.dat$variable=="y_m_o"|subpop.dat$variable=="o_m_o"|subpop.dat$variable=="y_f_o"|subpop.dat$variable=="o_f_o", "Other", "Hispanic" ))



diag.len <- nrow(diag.age.sex.dat)
msm.len <- nrow(msm.dat)
hiv.len <- nrow(subset(msm.dat[,c("pHIV_y","pHIV_o")], (!is.na(msm.dat[,"pHIV_y"])) & (!is.na(msm.dat[,"pHIV_o"]))))

plot.diag <- ggplot(diag.dat, aes(x=as.numeric(year), y=value, group=cat, color=cat)) + #plot reported early syphilis cases by age and sex
  geom_point() +
  geom_line( lwd=1) +
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Reported early syphilis cases\nper 100,000 population\n", title=statename) +
  coord_cartesian(ylim=c(0,75)) 

assign(paste("plot.diag", state, sep=""), plot.diag)

plot.pmsm <- ggplot(pmsm.dat, aes(x=as.numeric(year), y=value, color=cat)) + 
  geom_point(size=3) +
  geom_line(lwd=1) + 
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Proportion of male cases\nin MSM\n", title=statename) +
  coord_cartesian(ylim=c(0,1)) 
assign(paste("plot.pmsm", state, sep=""), plot.pmsm)

plot.hiv <- ggplot(hiv.dat, aes(x=as.numeric(year), y=value, color=cat)) + 
  geom_point(size=3) +
  geom_line(lwd=1) + 
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Proportion of MSM cases\nwith HIV coinfection\n") +
  coord_cartesian(ylim=c(0,1)) + 
  scale_x_continuous(breaks=seq((end.year-hiv.len+1),end.year,1))

assign(paste("plot.hiv", state, sep=""), plot.hiv)

plot.sec <- ggplot(sec.dat, aes(x=as.numeric(year), y=value, group=cat, color=cat)) + 
  geom_point() +
  geom_line( lwd=1) +
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Proportion of cases diagnosed\nwith secondary syphilis\n", title=statename) +
  coord_cartesian(ylim=c(0,1)) 
assign(paste("plot.sec", state, sep=""), plot.sec)

plot.el <- ggplot(el.dat, aes(x=as.numeric(year), y=value, group=cat, color=cat)) + 
  geom_point() +
  geom_line( lwd=1) +
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Proportion of cases diagnosed\nwith early latent syphilis\n") +
  coord_cartesian(ylim=c(0,1)) 

assign(paste("plot.el", state, sep=""), plot.el)


plot.subpop <- ggplot(subpop.dat, aes(x=as.numeric(year), y=value, color=subpop)) + 
  geom_point() +
  geom_line(lwd=1) +
  facet_wrap(sex~agecat) + 
  theme_classic()  +
  theme(legend.position="right",  axis.text.x=element_text(angle=90,hjust=0.5,size=12),axis.text.y=element_text(size=12),  
        plot.title=element_text(size=16, hjust=0),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  labs( x="Year", y="Reported early syphilis cases\nper 100,000 population\n") +
  coord_cartesian(ylim=c(0,180)) 
assign(paste("plot.subpop", state, sep=""), plot.subpop)
}

### save plots to pdf in Model_outputs folder ###
pdf(file="Model_outputs/syph_input_data_combined.pdf", width=10, height=8, paper="USr")
grid.arrange(plot.diagLA, plot.diagMA,
             plot.subpopLA, plot.subpopMA,
             ncol=2)
grid.arrange(plot.secLA, plot.secMA,
             plot.elLA,  plot.elMA,
             ncol=2)
grid.arrange(plot.pmsmLA, plot.pmsmMA,
             plot.hivLA,  plot.hivMA,
             ncol=2)
dev.off()


## plot diagnosed rates by age, sex, subpop ##
legend.diag <- get_legend(plot.diagLA)
plot.diagLA <- plot.diagLA + theme(legend.position ="none")
plot.diagMA <- plot.diagMA + theme(legend.position ="none")
legend.subpop <- get_legend(plot.subpopLA)
plot.subpopLA <- plot.subpopLA + theme(legend.position ="none")
plot.subpopMA <- plot.subpopMA + theme(legend.position ="none")

png(file="Model_outputs/syph_input_data_reported_rates.png", units="in", width=12, height=8, res=600)
grid.arrange(plot.diagLA, plot.diagMA, legend.diag,
             plot.subpopLA, plot.subpopMA, legend.subpop,
             ncol=9,nrow=3,
             layout_matrix = rbind(c(1,1,1,1,2,2,2,2,3), c(4,4,4,4,5,5,5,5,6), c(4,4,4,4,5,5,5,5,6)))
dev.off()

## plot proportion of cases in MSM and HIV+ MSM ##
legend.msm <- get_legend(plot.pmsmLA)
plot.pmsmLA <- plot.pmsmLA + theme(legend.position ="none")
plot.pmsmMA <- plot.pmsmMA + theme(legend.position ="none")
legend.hiv <- get_legend(plot.hivLA)
plot.hivLA <- plot.hivLA + theme(legend.position ="none")
plot.hivMA <- plot.hivMA + theme(legend.position ="none")

png(file="Model_outputs/syph_input_data_msm.png", units="in", width=12, height=8, res=600)
grid.arrange(plot.pmsmLA, plot.pmsmMA, legend.msm,
             plot.hivLA,  plot.hivMA,  blankPlot,
             ncol=9,nrow=2,
             layout_matrix = rbind(c(1,1,1,1,2,2,2,2,3), c(4,4,4,4,5,5,5,5,6)))
dev.off()
