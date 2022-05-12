library(faraway)

data <- gala

summary(lm(Species ~ Endemics + Elevation + Adjacent, data = data))

bayes_lm <- function(form, data, init, something = NULL, handle.na = NULL, 
                     n.iter = 1000) {
  
  m <- n.iter
  f <- form
  pb <- txtProgressBar(min = 1, max = m, title = "Process", style = 3)
  if (!is.list(init)) {stop("Initial values must be specified in a list.")}
  
  if      (is.null(handle.na))     {na <- "na.fail"}
  else if (handle.na == "na.omit") {na <- "na.omit"}
  
  d <- model.frame(f,data, na.action = na)
  x <- model.matrix(f,data = d)
  y <- matrix(model.frame(f,d)[,1], ncol = 1)
  n <- nrow(d)
  
  b0 <- 0; b1 <- 0; b2 <- 0; b3 <- 0; var <- 1
  
  #m.post <- lapply(1:m, FUN = function(i) {
  #  lapply(1:l)
  #})
  
  #for (i in 1:m) {
  #  m0.post <- ((sum(y - b[-1] %*% x[,-1])/var)+ mu/tau2)/((n/var)+1/tau2)
  #  b0[i+1] <- rnorm(1, m0.post, s0.post)
  #}
  
  s2 <- matrix(1,m)
  b  <- matrix(0,m,ncol(x))
  
  xtxi <- solve(t(x)%*%x)
  pars <- coef(lm(y~x-1))
  
  for (i in 2:m) {
    
    b[i,] <- pars+t(rnorm(ncol(x), mean=0, sd=1))%*%chol(s2[i-1]*xtxi)
    s2[i] <- 1/rgamma(1,nrow(x)/2, .5*t(y-x%*%(b[i,]))%*%(y-x%*%(b[i,])))
    
    setTxtProgressBar(pb, i)
    
  }
  cat("\n")
  
  coefs <- b
  re    <- s2
  
  return(list(coefficients = colMeans(b), residual.error = mean(s2), b, s2))
  
}

form <- as.formula(Species ~ Endemics + Elevation + Adjacent)
form2 <- as.formula(Species ~ Endemics + Nearest)

x <- model.matrix(form2, data)
y <- model.frame(form2, data)[,1]

n <- nrow(x)

lm.fit <- lm(y~x)

init <- c(0,0,1)
b0 <- -50
b1 <- 10
b2 <- 100
s2  <- 1

mu00 <- 0
tau00 <- 10000
mu10 <- 0
tau10 <- 10000
mu20 <- 0
tau20 <- 10000
A0 <- 0.0001
B0 <- 0.0001


out <- matrix(NA, nrow = 10000, ncol = 4)

for (i in 1:10000) {
  
  mu01 <- ((sum(y - (x[,2]*b1 + x[,3]*b2))/s2)+(mu00/tau00))/((n/s2)+(1/tau00))
  tau01 <- 1 / (n/s2+1/tau00)

  b0 <- rnorm(1, mu01, sqrt(tau01))

  mu11 <- ((sum(x[,2]*(y - (b0 + x[,3]*b2)))/s2)+(mu10/tau10)) / (sum(x[,2]^2)/s2+1/tau10)
  tau11 <- 1 / (sum(x[,2]^2)/s2 + 1 / tau10)

  b1 <- rnorm(1, mu11, sqrt(tau11))
  
  mu21 <- ((sum(x[,3]*(y - (b0 + x[,2]*b1)))/s2)+(mu20/tau20)) / (sum(x[,3]^2)/s2+1/tau20)
  tau21 <- 1 / (sum(x[,3]^2)/s2 + 1 / tau20)
  
  b2 <- rnorm(1, mu21, sqrt(tau21))
  
  A1 <- n/2 + A0
  B1 <- sum((y - (b0 + b1*x[,2] + b2*x[,3]))^2)/2 + B0
  s2  <- 1 / rgamma(1, A1, B1)
  
  out[i,] <- c(b0,b1,b2,s2)
  
}

round(colMeans(out[-(1:100),]), 5)

plot(out[-c(1,2,3),4], type = "l")

sqrt(869.99887)

sqrt(840.20)

869.17 * (29/30)

out[3,]

plot(out[,4], type = "l")

head(out, 20)

all.equal(head(out[,4], 100), out[101:200,4])

var(out)

dat <- MASS::birthwt
summary(lm(bwt ~ age + lwt, data = dat))
var(out)

lm.fit$fitted.values

summary(lm.fit)


mean(1 / rgamma(10000, n/2, .5*t(y - x%*%c(b0,b1,b2))%*%(y - x%*%c(b0,b1,b2))))
28^2
plot(density(x = 1/rgamma(1000, .1, .1)))
lines(density(x = 1/rgamma(1000, 1,1)), col = "blue")
lines(density(x = 1/rgamma(1000, 10,10)), col = "red")
lines(density(x = 1/rgamma(1000, 100,100)), col = "purple")
lines(density(x = 1/rgamma(1000, 100000,100000)), col = "purple")
lines(density(x = 1/rgamma(1000, 1/100,1/100)), col = "purple")

hist(rgamma(1000, .01, .01))
plot(density(rgamma(1000, 10000, 10000)))


rnorm(1000, 0, sd =10000000)


out[1:100,]

A1 <- n/2 + A0
B1 <- sum((y - (b0 + b1*x[,2]))^2)/2 + B0


b <- init[1:2]
s <- init[3]

numsim <- 100000
samples <- matrix(0, numsim, 3)

for (i in 1:numsim) {
  VVV <- n/s + 1/.001
  MMM <- sum(y - x[,2]*b[2])
  b[1] <- rnorm(1, MMM/VVV, 1/sqrt(VVV))

  VVV <- sum(x[,2]^2)/s + 1/.001
  MMM <- sum(x[,2]*(y-b[1]))/s
  b[2] <- rnorm(1, MMM/VVV, 1/sqrt(VVV))
  
  SSE <- sum((y-b[1]-x[,2]*b[2])^2)
  s <- 1 / rgamma(1, n/2 + .001, SSE/2+ .001)
  
  samples[i,] <- c(b,s)
}
colMeans(samples)
head(samples)
out <- Gibbs.regression(x[,2], y, M = NULL, Nsamples = 1000, trace = 'bsmt', fix = 'xy', intercept = TRUE)

plot(1:100,
     rgamma(100, shape = 10000, scale = 10000))

Gibbs.regression

out <- bayes_lm(form, data = data, init = a, n.iter = 100000)
out[[1]]

out

x1 <- rnorm(100)
x2 <- x1 + rnorm(100,0,.2)
cor(x1,x2)

y1 <- x1 + x2
y2 <- x1-x2

cor(y1,y2)

x <- c(1,2,3)
names_x <- NA

for (i in 1:length(x)) {
  names_x[i] <- paste("beta",i-1, sep = "")
}
names(x) <- names_x
x

names(x) <- sapply(x, FUN = function(y) paste("beta", y, sep = ""))
