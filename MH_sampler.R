x_times_2_sampler <- function(initial) {
  if (abs(initial) > 1) stop("The initial value must lay between 0 and 1")
  x <- runif(1)
  r <- ((2*x) / (2*initial))
  if (r >= runif(1)) return(x)
  else return(initial)
}

out <- numeric(length = 100000)

x <- .5

for (i in 1:100000) {
  x <- x_times_2_sampler(initial = x)
  out[i] <- x
}

hist(out,freq=F)
