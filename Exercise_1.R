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

plot(density(res[[2]][,1]))

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
## Theta_PE = (1 + 40)/(1+40+1+143-40) = (41/145) = .2828
## RR = .4126 / .2828 = .68533

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
## Since I do not think the data is as informative as it could be, due
## to the fact that there are large differences between all studies, 
## I think it would be wise to use the previous studies as power priors,
## both with a prior weight of .05. 
## Then, the alpha would become .05 * (120 + 40)  + 1 = 9 for PE,
## the beta would become .05 * ((245+105)-(120+40)) + 1 = 10.5 for PE.
## And the alpha for PC would become .05 * (80+45) + 1 = 7.25,
## and the beta for PC would become .05 * ((275+130)-(80+45)) + 1 = 15.

## Let's now first plot the four beta functions, to get a general idea
## how the data is distributed. The solid lines indicate the PC condition,
## while the dashed lines indicate the PE condition. Furthermore, pink
## and purple indicate the female veteran study, while the blue lines
## indicate the male veteran study.
p <- seq(0,1, length = 1000)
plot(p, dbeta(p, 120, 125), type = "l", lty = 2, col = "pink")
lines(p, dbeta(p, 80, 195), type = "l", col = "purple")
lines(p, dbeta(p, 40, 65), type = "l", lty = 2, col = "blue")
lines(p, dbeta(p, 45, 85), type = "l", col = "light blue")

model.def_J <- jags.model(file = 'Exercise 1 - Model_J.txt', data = dat, n.chains = 2)

update(object = model.def_J, n.iter = 5000)

par_J <- c('theta.PE', 'theta.PC', 'RR')
res_J <- coda.samples(model = model.def_J, variable.names = par, n.iter = 500000)

summary(res)
summary(res_J)

## The results are barely influenced, which is of course predominantly 
## due to the fact that I did not regard the prior studies as very 
## relevant. However, mainly because of the learning potential, I wanted
## to work with both studies. What can be seen though, is that the 
## relative risk slightly decreased, although the difference occurs in
## the third decimal. However, in substantive terms, this means that
## the relative difference between PC and PE has become even larger, 
## indicating that people receiving the PC treatment have an estimated
## probability of losing the PTSD diagnosis that is .69 times the
## probability of losing the PTSD diagnosis for patients receiving the
## PE treatment. Furthermore, the 95% credibility interval becomes 
## somewhat narrower, from [.4891-.9422] to [.5007-.9221], indicating
## that with 95% certainty, the value of the parameter given the data 
## and given the prior lies somewhat further away from zero, accentuating
## the better prognosis after PE.

## K
model.def_K <- jags.model(file = 'Exercise 1 - Model_K.txt', data = dat, n.chains = 2)

update(object = model.def_K, n.iter = 5000)

par_K <- c('theta.PE', 'theta.PC', 'RR')
res_K <- coda.samples(model = model.def_K, variable.names = par, n.iter = 500000)

summary(res)
summary(res_J)
summary(res_K)

## When we give both prior studies a weight of .5 instead of .05,
## we find that the mean of the RR remains more or less unchanged (it 
## decreases with approximately .008 over the three models). However,
## the mean difference is the narrowing of the CI, from [.4891-.9422]
## in the first model to [.5532-.8334] in the last model, indicating
## that the precision regarding the estimate of the parameter value
## increases.

## L

## theta_pc_j = .2859
## theta_pe_j = .4174
## RR_j = .6850
## theta_pc_k = .2978
## theta_pe_k = .4371
## RR_k = .6813

## The analytically derived results are highly similar to the 
## simulated results.