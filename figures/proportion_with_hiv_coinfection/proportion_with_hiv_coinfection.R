
devtools::load_all()
library(ggplot2)

state <- 'LA'

load.start()

la_df <- data.frame(state = 'Louisiana', age = rep(c('20-44 y', '45-64 y'), each=3), prop = p.hiv.cases)

state <- 'MA'

load.start()

ma_df <- data.frame(state = 'Massachusetts', age = rep(c('20-44 y', '45-64 y'), each=5), prop = p.hiv.cases)


df <- rbind.data.frame(la_df, ma_df)
df$year <- c(rep(2014:2016, 2), rep(2012:2016, 2))

ggplot(df, aes(x = year, y = prop*100, color = age, shape = age)) + 
  geom_point(size = 5, alpha = 0.8) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(~state) + 
  expand_limits(y=c(0,100)) + 
  ylab("Percentage Among MSM with HIV Coinfection (%)") + 
  xlab("Years") + 
  ggtitle("Percentage of Cases Among MSM with HIV Coinfection") + 
  scale_color_manual(values=c('brown2', 'cadetblue')) + 
  scale_shape_manual(values=c(16,15)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,vjust=.65)) 


ggsave("proportion_with_hiv_coinfection.png", height = 4.5, width = 8)
