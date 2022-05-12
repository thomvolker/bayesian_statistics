model{
  for(i in 1:length(watts)) {
    watts[i] ~ dnorm(mu,tau)
    mu[i] <- b[1] + b[2]*speed[i] + b[3]*distance + b[4]*duration
  }
  for(k in 1:4) {
    b[k] ~ dnorm(0,.001)
  }
  tau ~ dgamma(.001,.001)
  sigma <- 1/tau
}