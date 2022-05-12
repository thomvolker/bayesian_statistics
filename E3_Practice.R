source("Exercise 3 - Data.txt")

md <- jags.model("E3_Practice.txt", data = dat, n.chains = 2)

update(md,n.iter = 100)

par <- c("theta.PE", "theta.PC")

samples <- coda.samples(md, par, n.iter=10000)

plot(samples)

summary(samples)

proportion.half1.PE <- mean(dat$LD.PE[1:70])
proportion.half2.PE <- mean(dat$LD.PE[71:141])

proportion.half1.PC <- mean(dat$LD.PC[1:71])
proportion.half2.PC <- mean(dat$LD.PC[72:143])

theta.PE.chain1 <- samples[[1]][,"theta.PE"]
theta.PE.chain2 <- samples[[2]][,"theta.PE"]
theta.PC.chain1 <- samples[[1]][,"theta.PC"]
theta.PC.chain2 <- samples[[2]][,"theta.PC"]

replicates <- list(NA)

replicates[[1]] <- array(NA, dim = c(length(theta.PE.chain1), dat$n.PE))
replicates[[2]] <- array(NA, dim = c(length(theta.PE.chain2), dat$n.PE))
replicates[[3]] <- array(NA, dim = c(length(theta.PC.chain1), dat$n.PC))
replicates[[4]] <- array(NA, dim = c(length(theta.PC.chain2), dat$n.PC))

pb <- txtProgressBar(min = 0, max = length(theta.PE.chain1), char = "+")


for(t in 1:length(theta.PE.chain1)) {
  replicates[[1]][t,] <- rbinom(dat$n.PE,1,theta.PE.chain1[t])
  replicates[[2]][t,] <- rbinom(dat$n.PE,1,theta.PE.chain2[t])
  replicates[[3]][t,] <- rbinom(dat$n.PC,1,theta.PC.chain1[t])
  replicates[[4]][t,] <- rbinom(dat$n.PC,1,theta.PC.chain2[t])
  setTxtProgressBar(pb,t)
}

replicates[[1]]

theta.PE.dif <- proportion.half1.PE - proportion.half2.PE
theta.PC.dif <- proportion.half1.PC - proportion.half2.PC

difs <- c(theta.PE.dif,theta.PE.dif,
          theta.PC.dif,theta.PC.dif)


dif_means <- function(x) {
  dif <- rowMeans(x[,1:floor(ncol(x)/2)]) - rowMeans(x[,ceiling(ncol(x)/2):ncol(x)])
  return(abs(dif))
}
mean(difs[1] > dif_means(replicates[[2]]))
out <- NA
for (i in 1:4) out[i] <- mean(lapply(replicates, FUN = function(x) dif_means(x))[[i]] > abs(difs)[i])
out
lapply(replicates, FUN = function(x) colMeans(x))

r <- runif(1)

eq <- replicate(100000000, expr = r == runif(1))
sum(eq)
head(eq)

curve(dnorm(x, 4, 2),xlim = c(-4,12))
curve(dnorm(x, 2, 2),add = TRUE, col = "blue")
abline(v = 2)
abline(v = 4, col = "blue")
abline(h = dnorm(2,4,2))
