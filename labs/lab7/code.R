# data from http://johnstonsarchive.net/other/worldpop.html
pop <- read.csv("data/pop2.csv")

japan <- pop[,c("year","Japan")] %>% mutate(Japan = as.numeric(Japan)) %>% 
  subset(!is.na(Japan)) %>% subset(year < 1901 | year > 1945)
plot(japan)


# logistic model

N.logistic <- function(x, N0) K/(1 + ((K - N0)/N0)*exp(-r0*x))



