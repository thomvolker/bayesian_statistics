source("Exercise 4 - Data.txt")

dat <- as.data.frame(dat)

head(dat)

library(rjags)

model.def <- jags.model(file = "Exercise 4 - Model.txt", data = dat, n.chains = 2)

update(model.def, n.iter = 1000)

summary(glm(depr ~ 1 + gender, data = dat, family = binomial))

par <- c("a", "b", "OR")

res <- coda.samples(model = model.def, variable.names = par, n.iter = 10000)

summary(res)

traceplot(res)
autocorr.plot(res)

densplot(res)
gelman.plot(res)

## Convergence immediately established.

## The estimated odds ratio
HPDinterval(mcmc(do.call(rbind, res)))
## The 95% highest posterior density interval ranges from
## 1.07 to 2.33, so the posterior probability that females
## have a higher predicted probability to develop depression
## is at least 95%. Also, the predicted odds of females to
## develop depression are on average 1.66 times the odds of
## males. And in terms of probabilities, the results 
## indicate the following:
exp(sum(summary(res)[[1]][2:3,1])) / 
  (1 + exp(sum(summary(res)[[1]][2:3,1])))
## .17 for females, and
exp(summary(res)[[1]][2,1]) / 
  (1 + exp(summary(res)[[1]][2,1]))
## .11 for males.

model.def.E <- jags.model(file = "Exercise 4_E - Model.txt",
                          data = dat, n.chains = 2)

update(model.def.E, n.iter = 1000)

res.E <- coda.samples(model = model.def.E, variable.names = par, n.iter = 10000)

summary(res.E)
traceplot(res.E)
gelman.plot(res.E)

dat$depr <- as.factor(dat$depr)


summary(glm(depr ~ 1 + as.factor(gender) + as.factor(cidi1) + as.factor(depGP) + as.factor(antidep), data = dat,
        family = "quasibinomial"))

DIC <- dic.samples(model.def, n.iter = 10000)
DIC.E <- dic.samples(model.def.E, n.iter = 10000)
dic.sa
DIC


x1 <- rnorm(100, 5, 7)
x2 <- x1 + rnorm(100, 0, 3)

cor(matrix(c(x1,x2),nrow = 100,ncol=2))

var(matrix(c(x1,x2),nrow = 100,ncol=2))

x1c <- (x1 - mean(x1))/sd(x1)
x2c <- (x2 - mean(x2))/sd(x2)

var(matrix(c(x1c,x2c),nrow = 100,ncol=2))

plot(rgamma(100, .1,.1))


quantile(rbinom(1000000, 422, prob = .5), c(.05,.95))

combn(1:4, 3)

var(rbinom(1000, 1000, .5))

normal <- function(x, mu, sigma) {
  1 / sqrt(2*pi*sigma^2) * exp(-((x-mu)^2 / (2*sigma^2)))
}

curve(normal(x,10,5), -20,40)

normal(-10:30, 10,5)

x0 <- rep(1,100)
x1 <- rnorm(100, 2)
x2 <- rnorm(100, 3)
y <- x1 + x2 + rnorm(100)

summary(lm(y ~ x1 + x2))

b <- solve(t(cbind(x0,x1,x2))%*%cbind(x0,x1,x2))%*%(t(cbind(x0,x1,x2))%*%y)
s <- 1/100 * t(y - cbind(x0,x1,x2)%*%b)%*%(y - cbind(x0,x1,x2)%*%b)

acov <- as.numeric(s)*(solve(t(cbind(x0,x1,x2))%*%cbind(x0,x1,x2)))
sqrt(diag(acov))

?solve
x_2 <- function (x) 2*x
x <- seq(-10,30, by = .00001)
curve(x_2(x))
sample(x, prob = normal(x, 10, 5))