dat <- source("Exercise 4 - Data.txt")

dat <- as.data.frame(dat$value)

init1 <- list(b = c(-10,10))
init2 <- list(b = c(5,15))
inits <- list(init1, init2)

md <- jags.model("E4_model.txt", data = dat, n.chains = 2, inits = inits)

update(md, n.iter = 1000)

par <- c("b", "OR")

samples <- coda.samples(md, par, n.iter = 1000)
dic.samples(md, n.iter=5000)
plot(samples)
autocorr.plot(samples)
gelman.plot(samples)
densplot(samples)

init1 <- list(b = c(-3,-5,5,8,3))
init2 <- list(b = c(5,8,5,-3,5))
inits <- list(init1, init2)

md_E1 <- jags.model("E4_model_E1.txt", data = dat, n.chains = 2)

update(md_E1, n.iter = 5000)

par <- c("b", "OR")

dic1 <- dic.samples(md_E1,n.iter = 10000)
samples <- coda.samples(md_E1, par, n.iter = 1000)
plot(samples)

md_E2 <- jags.model("E4_model_E2.txt", data = dat, n.chains = 2)

update(md_E2, n.iter = 5000)

par <- c("b", "OR")

dic2 <- dic.samples(md_E2,n.iter = 5000)
samples <- coda.samples(md_E2, par, n.iter = 1000)
plot(samples)
summary(samples)
HPDinterval(samples)

x <- matrix(1:9,3,3)
y <- c(1,2,0)
x%*%y

x%*%x^(-1)
solve(x)
t(diag(4))

sum(solve(matrix(c(1,2,2,1),2,2)))
sum(solve(matrix(c(1,.02,.02,1),2,2)))

curve(dbeta(x,2,1))
