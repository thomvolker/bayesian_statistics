library(foreign)
library(rjags)

dat <- read.spss("Exercise 2 - Data.sav", to.data.frame = TRUE)

head(dat)

init1 <- list(b = c(0,0,0),
              tau = 1)
init2 <- list(b = c(-10,10,5),
              tau = 10)
inits <- list(init1, init2)

md <- jags.model("E2_Practice.txt", data = dat, inits = inits, n.chains = 2)

update(md, n.iter = 500)

par <- c("b", "sigma")

res_1 <- coda.samples(md,par,n.iter = 1000)

summary(res)

traceplot(res)
autocorr.plot(res)
densplot(res)
gelman.plot(res)

plot(res)

init1_2 <- list(b = c(0,0,0,0),
              tau = 1)
init2_2 <- list(b = c(-10,10,5,20),
              tau = 10)
inits_2 <- list(init1_2, init2_2)

md <- jags.model("E2_Practice_2.txt", data = dat, inits = inits_2, n.chains = 2)

update(md, n.iter = 2000)

par <- c("b", "sigma")

res <- coda.samples(md,par,n.iter = 20000)

summary(res)

plot(res)

res_1
