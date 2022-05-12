fun <- function(x) 2*x

last <- 1
val <- NA

numsim <- 2000

for (i in 1:numsim) {
  candidate <- runif(1)
  alpha <- fun(candidate)/fun(last)
  if (runif(1) < min(alpha, 1)) last <- candidate
  val[i] <- last
}

hist(val)

plot(val, type = "l")

target <- function(x) {
  return(ifelse(x<0,0,exp(-x)))
}

target(10)
exp(10)

x<-10

for(i in 2:100000) {
  current_x <- x[i-1]
  proposed_x <- current_x + rnorm(1)
  A <- target(proposed_x) / target(current_x)
  if (runif(1) < A) x[i] <- proposed_x
  else x[i] <- current_x
}

plot(x, type = "l")

hist(x)
