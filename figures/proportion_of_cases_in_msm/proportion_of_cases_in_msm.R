devtools::load_all()
library(ggplot2)

state <- 'LA'

load.start()

la_df <- data.frame(state = 'Louisiana', age = rep(c('20-44 y', '45-64 y'), each=5), prop = p.msm.cases)

state <- 'MA'

load.start()

ma_df <- data.frame(state = 'Massachusetts', age = rep(c('20-44 y', '45-64 y'), each=5), prop = p.msm.cases)


df <- rbind.data.frame(la_df, ma_df)
df$year <- rep(2012:2016, 4)

ggplot(df, aes(x = year, y = prop*100, color = age, shape = age)) + 
  geom_point(size = 5, alpha = 0.8) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(~state) + 
  expand_limits(y=c(0,100)) + 
  ylab("Percentage Among MSM (%)") + 
  xlab("Years") + 
  ggtitle("Percentage of Male Cases Among MSM") + 
  scale_color_manual(values=c('orange', 'royalblue')) + 
  scale_shape_manual(values=c(16,15)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,vjust=.65)) 


ggsave("proportion_of_cases_in_msm.png", height = 4.5, width = 8)
