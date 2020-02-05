library(rjags)

## A 
source('Exercise 1 - Data.txt')

dat

## B
model.def <- jags.model(file = 'Exercise 1 - Model.txt', data = dat, n.chains = 2)

## C
## Not neccessary

## D
update(object = model.def, n.iter = 5000)

par <- c('theta.PE', 'theta.PC', 'RR')
res <- coda.samples(model = model.def, variable.names = par, n.iter = 500000)

## E
## Not necessary

## F

summary(res)

quantile(res[[2]][,1], probs = seq(.9,1, .01))

## The result is quite convincing in favor of PE, due to the fact that the 
## relative risk of the loss of diagnosis after Prolonged Exposure is .68
## as compared to Present-Centered therapy. Furthermore, the probability of
## losing the diagnosis of PTSD for people that receive Prolonged Exposure
## is approximately .41, while the probability of losing the PTSD diagnosis
## is .28 for people who receive Present-Centered therapy. Also, the 95%
## credibility interval of the OR is [.49 - .94], and since 1 is not included
## in this interval, the probability that the effect PC better than or equal to
## the effect of PE is smaller than 1%.

## G

## Theta_PC = (1 + 58)/(1+58+1+141-58) = (59/143) = .4126
## Theta_PE = (1 + 40)/(1+40+1+141-40)

## H

## I would say that there is some prior information in the other datasets,
## that is, those could be used as a primer for the current study. However,
## table 2 concerns females, who may react to the different therapies in 
## different ways due to biological or behavioral differences between
## genders. Furthermore, no date is specified for this data, and the
## corresponding research procedure isn't either. Therefore, if we assume
## that the research was done properly and not too long ago, we might use
## it as a power prior with a weight of .05. The second study, with the
## WWII data is very long ago, and therefore not extremely informative,
## since it might well be that the procedure has changed since then. In
## favor of this data is however, that it concerns all male veterans,
## which also holds for the current study. Therefore, in my opinion we
## could use this study as a power prior as well, also with a weight 
## of .05.

## I
