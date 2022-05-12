#########################
## Question A
#########################

## The binomial distribution assumes that every person within a
## condition has the same probability of losing the diagnosis which
## might not be the case, since this probability might depend on
## other factors, such as severity of the disease, stadium, and
## maybe other demographical characteristics.

#########################
## Question B + C
#########################
library(rjags)

source("Exercise 3 - Data.txt")

dat

init1 <- list(theta.PE = 1, theta.PC = 1)
init2 <- list(theta.PE = 0, theta.PC = 0)
init <- list(init1, init2)

model.def <- jags.model(file = 'Exercise 3 - Model_template.txt', inits = init, data = dat, n.chains = 2)

update(object = model.def, n.iter = 1000)

par <- (c("theta.PE", "theta.PC"))

res <- coda.samples(model = model.def, variable.names = par, n.iter = 100000)

#summary(res)

#plot(res)

#########################
## Question D
#########################

## Differences in the two halves might indicate that there
## is a systematic difference between the first sampled
## persons and persons who are sampled later. If this is 
## the case, both groups should have different probabilities
## and they aren't independently distributed either.

# test statistic for the PE condition 
proportion.half1.PE <- sum(dat$LD.PE[1:70])/70 
proportion.half2.PE <- sum(dat$LD.PE[71:141])/71 
diff.PE <- proportion.half1.PE - proportion.half2.PE 
# test statistic for the PC condition 
proportion.half1.PC <- sum(dat$LD.PC[1:71])/71 
proportion.half2.PC <- sum(dat$LD.PC[72:143])/72 
diff.PC <- proportion.half1.PC - proportion.half2.PC

#########################
## Question E
#########################

theta.PE.chain1 <- res[[1]][,"theta.PE"] 
theta.PE.chain2 <- res[[2]][,"theta.PE"] 
theta.PC.chain1 <- res[[1]][,"theta.PC"] 
theta.PC.chain2 <- res[[2]][,"theta.PC"]
theta.PE.chain1

library(LaplacesDemon)

replicates <- list(NA)

replicates[[1]] <- array(data = NA, dim = c(length(theta.PE.chain1), dat$n.PE)) 
replicates[[2]] <- array(data = NA, dim = c(length(theta.PE.chain2), dat$n.PE)) 
replicates[[3]] <- array(data = NA, dim = c(length(theta.PC.chain1), dat$n.PC)) 
replicates[[4]] <- array(data = NA, dim = c(length(theta.PC.chain2), dat$n.PC))



for(t in 1:length(theta.PE.chain1)) {
  # ... sample a replicated dataset by sampling n times from the Bernoulli distribution 
  # with as probability of success on each trial the parameter estimate 
  replicates[[1]][t,] <- rbern(n = dat$n.PE, prob = theta.PE.chain1[t]) 
  replicates[[2]][t,] <- rbern(n = dat$n.PE, prob = theta.PE.chain2[t]) 
  replicates[[3]][t,] <- rbern(n = dat$n.PC, prob = theta.PC.chain1[t]) 
  replicates[[4]][t,] <- rbern(n = dat$n.PC, prob = theta.PC.chain2[t])
}

dim(replicated.PC.chain1)
lapply(replicates, FUN = function(x) rowMeans(x))

dif_means <- function(rep1) {
  
  rep_half1 <- rowMeans(rep1[,1:round(ncol(rep1)/2)])
  rep_half2 <- rowMeans(rep1[,round(ncol(rep1)/2):ncol(rep1)])
  
  mean_dif <- rep_half1 - rep_half2
  
  return(abs(mean_dif))
}


out <- NA
for (i in 1:4) out[i] <- mean(lapply(replicates, FUN = function(x) dif_means(x))[[i]] > abs(c(diff.PE, diff.PE, -diff.PC, -diff.PC)[i]))
out


