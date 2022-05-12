library(foreign)
library(rjags)

dat <- read.spss("Exercise 2 - Data.sav", to.data.frame = TRUE)

head(dat)

## jags.model might need initial values, in case of multiple chains, so will also
## have to specify initial values by means of the "init" command

init1 <- list(alpha = 10, beta1 = 100, beta2 = -100, beta3 = -10, tau = 5)
init2 <- list(alpha = 100, beta1 = -100, beta2 = -10, beta3 = 10, tau = .0001)
init <- list(init1, init2)

list(0,1,1)

model.def <- jags.model(file = 'Exercise 2 - Model_Int.txt', 
                        inits = init, data = dat, n.chains = 2)

update(object = model.def, n.iter = 1)

par <- (c("alpha", "beta1", "beta2", "beta3", "tau"))

res <- coda.samples(model = model.def, variable.names = par, n.iter = 100000)

plot(res)

b0 <- summary(res)[[1]][2,1]
b1 <- summary(res)[[1]][3,1]
b2 <- summary(res)[[1]][4,1]

traceplot(res)
autocorr.plot(res)
densplot(res)
gelman.plot(res)

## also look at the monte carlo error (naive SE)



## The intercept is estimated to be .868, indicating that persons scoring
## 0 on the variables extraversion and agreeableness, have a predicted value
## of .868 on attitude towards pet ownership.

sum(res[[1]][,2] > 0)/50000

## The regression coefficient of extraversion is estimated to be -.082, 
## indicating that a 1 point increase on extraversion is associated with
## an -.082 point decrease on attitude towards pets. Furthermore, the 95%
## credibility interval ranges from -.1936 to .02, indicating that there 
## is a 95% probability that the true population value is located in this
## interval. Then, 7.3 percent of the sampled values are bigger than 0,
## indicating that there is a 7.3 percent probability that the true value
## of this regression coefficient is bigger than 0 (although this would)
## be a minor effect.

sum(res[[1]][,3] > 0)/50000

## The regression coefficient of agreeableness is estimated to be .26,
## indicating that a 1 point increase in agreeableness is related to a
## .26 point increase in attitude towards pet ownership. Furthermore, the
## 95% credibility interval ranges from .08 to .43, while the proportion 
## of sampled values that is bigger than 0 is .998, indicating that there
## is a 99.8% probability that the true regression coefficient for the 
## variable agreeableness is bigger than zero.

sum()

## Than, lastly, we have to look at the interaction in the model.
summary(res)
## The interaction regression coefficient is estimated to be .11,
## indicating that a higher score on both explanatory variables is 
## related to a higher score on attitude towards pets. However, when
## one of the two is negative, this effect will become negative, and 
## would thus lead to a decrease in attitude towards pets. If both IV's 
## are negative, this effect is probably also negative.


dat$extraversion


autocorr.plot(res)


