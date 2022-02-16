
rm(list=ls())
WA <- read.delim("/mnt/box/Box Sync/teaching/EFB370 - Population Ecology/data/WA_SeaOtters_PopGrowth.dat")
plot(WA, log = "y")
head(WA)

WA.growth <- lm(log(count)~year, data = WA)
summary(WA.growth)
plot(WA.growth)
