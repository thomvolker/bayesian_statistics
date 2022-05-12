## Exercise 1
library(rjags)
source("Exercise 1 - Data.txt")

md <- jags.model("E1_Practice.txt", data = dat, n.chains = 2)

update(object = md, n.iter = 50)

par <- c("theta.PE", "theta.PC", "RR")

res <- coda.samples(md, par, 1000)

traceplot(res)
autocorr.plot(res)

summary(res)

md_2 <- jags.model("E1_Practice_2.txt", data = dat, n.chains = 2)

update(object = md_2, n.iter = 50)

par <- c("theta.PE", "theta.PC", "RR")

res <- coda.samples(md_2, par, 1000)

traceplot(res)
autocorr.plot(res)

summary(res)
