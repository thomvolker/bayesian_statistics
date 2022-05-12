library(rjags)
library(R2jags)

params <- c("b", "sigma")

fit_r2jags <- jags(na.omit(dat), parameters.to.save = params,
                   model.file = "jags_model_assignment.txt",
                   n.chains=6,n.iter=250,n.burnin=2)

bayes_lm(form = watts ~ nact + speed + distance, data = dat, prior = "uninf",
         DIC.samples = F, 
         init = list(init1 = list(b0 = 0, b1 = 0, b2 = 0, b3 = 0, sigma_squared = 1),
                     init2 = list(b0 = 100, b1 = 100, b2 = 100, b3 = 100, sigma_squared = 100)),
         burnin = 5, n.iter = 500, n.chains = 2, 
         seed = 12345, trace.plot = F, autocor.plot = F)
fit_r2jags
dat2 <- na.omit(dat[,-1])

mf <- jags.model("jags_model_assignment.txt", data = dat2, n.chains = 3)

model_txt <- "
model{
  for(i in 1:length(score)){
    score[i] ~ dnorm(mu1[i],tau)
    mu1[i] <- b0 + b1*group[i]
  }
  b0 ~ dnorm(0,.0001)
  b1 ~ dnorm(0,.00001)
  tau ~ dgamma(.0001,.0001)
  sigma <- 1/tau
  }"

practice_dat2
write.table(model_txt, row.names = F, col.names = F)
mf2 <- jags.model(textConnection(model_txt), data = practice_dat2, n.chains = 3) 

update(mf, n.iter = 1000)
update(mf2, n.iter = 1000)

params.to.save <- c("b", "sigma")
params.to.save2 <- c("b0", "b1", "sigma")

samples <- coda.samples(mf, params.to.save, n.iter = 10000)
samples2 <- coda.samples(mf2, params.to.save2, n.iter = 1000)
summary(samples)
HPDinterval(samples2)
class(samples)
gelman.diag(samples)
gelman.plot(samples)
plot(samples)
summary(samples2)
gd <- gelman.diag(samples)
gelman.diag(samples2)
gelman.plot(samples2)
gd$mpsrf
class(summary(samples))

summary(lm(data$score ~ 1))

summary(samples)

autocorr.plot(samples)

plot(samples)
summary(fit)


dic.samples(mf2, n.iter = 10000)
thesis_out$DIC
apply(thesis_out$posterior,2,sd)

