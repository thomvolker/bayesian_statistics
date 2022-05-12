############################################################################
## First part, data thingies
############################################################################
## Read in data (this is the full data, later, we will subset the data)
## and only use the first 88 rows.
data <- read.csv(file = "/Users/thomvolker/Downloads/export_34188992(1)/activities.csv")

# Load tidyverse for data tidying, lubridate for working with dates, and RColorBrewer for plots
library(tidyverse); library(lubridate); library(RColorBrewer); library(e1071)

dat_uncentered <- data %>% 
  ## First, we make the dates in the data as POSIXct objects
  mutate(date = mdy_hms(as.character(.$Activity.Date))) %>%
  ## Only keep activities that are "bikerides"
  filter(Activity.Type == "Ride") %>% 
  ## select variables of interest (date, duration, distance, watts)
  select(date = date, duration = Elapsed.Time.1, distance = Distance.1, watts = Average.Watts) %>%
  ## transform duration from seconds to hours, and distance to km and compute speed
  mutate(duration = duration / 3600, distance = distance / 1000, speed = distance / duration) %>%
  ## filter some erroneously included activities
  filter(distance != 0 | speed != 0)

## Also, I want to compute the number of activities in the four weeks 
##(28 days) before the activity of interest
dat_uncentered$nact <- NA
for (i in 1:nrow(dat_uncentered)) {
  from <- dat_uncentered$date[i] - (60*60*24*28) ## the date minus four weeks
  to   <- dat_uncentered$date[i]                 ## until the date of interest
                                      ## and sum these dates inbetween
  dat_uncentered$nact[i] <- sum(dat_uncentered$date > from & dat_uncentered$date < to, na.rm = TRUE)
}

## Center all variables (except the date variable, which will not be used anyway)
dat <- dat_uncentered[1:88,]
dat[,c(2,3,5,6)] <- as.data.frame(lapply(dat[,c(2,3,5,6)], FUN = function(x) x-mean(x, na.rm = TRUE)))

############################################################################
## Second part - Gibbs sampling
############################################################################

bayes_lm <- function(form, data, init = NULL, prior = "uninf", n.chains,  burnin, 
                     which.MH = NULL, MH.inital = NULL, MH.conditional = NULL, MH.proposal = NULL, 
                     handle.na = NULL, n.iter = 1000, conv.check = FALSE, DIC.samples = FALSE,
                     trace.plot = FALSE, autocor.plot = FALSE, DIC = FALSE, seed = NULL) {
  
  call_object <- match.call()
  
  if (n.chains > 9) {
    warning("Maximum number of chains is nine, so sampling continued with nine chains.")
    n.chains <- 9
  }
  
  set.seed(seed) #set seed if required for reproducibility

  x <- model.matrix(form, data) #matrix with predictors
  y <- model.frame(form, data)[,1] #column vector with the outcome
  
  p <- ncol(x) + 1; n <- nrow(x); #specify #parameters (including the residual variance) and sample size
  a_names <- list(NULL, c(colnames(x), "Residual variance"), #names for the arrays containing results
                     sapply(1:n.chains, function(x) paste("Chain", x)))
  
  pe_bi <- array(NA, dim = c(burnin, p, n.chains), dimnames = a_names) #empty matrix for burn in
  pe_ni <- array(NA, dim = c(n.iter, p, n.chains), dimnames = a_names) #empty matrix for sampled estimates
  
  pb_bi <- txtProgressBar(min = 1, max = burnin, char = "~", style = 3) #progress bar for burn in 
  pb_ni <- txtProgressBar(min = 1, max = n.iter, char = "~", style = 3) #progress bar for sampling
  pb_cp <- pb_ap <- txtProgressBar(min = 1, max = p,  char = "~", style = 3) #progress bar for trace plot and autocor plot
  
  if (length(prior) == 1) {if (prior == "uninf") { # if uninformative prior is used, specify
    m00 <- rep(0, p-1); t2_00 <- rep(10000, p-1)   # uninformative parameters inside the function
    A00 <- .0001; B00 <- .0001
  }}
  else if (length(prior) != 3) stop("Prior distribution should be specified in three vectors (mean, standard deviations, variance parameters).")
  else {
    if (length(prior[[1]]) != p-1 | length(prior[[2]]) != p-1 | length(prior[[3]]) != 2) {
      stop("Prior means and standard deviations should be specified for all regression coefficients including the intercept,
            and the variance parameters should be specified for both the shape and rate parameter.")
    }
    else { # if informative prior is used, this format should be used
      m00 <- prior[[1]]; t2_00 <- prior[[2]] # means and sd's
      A00 <- prior[[3]][1]; B00 <- prior[[3]][2] # and parameters for the variance
    }
  }
  
  if (is.null(init)) { #if no initial values are specified, draw initial values automatically
    init <- lapply(1:n.chains, FUN = function(x) c(rnorm(p-1,mean=0, sd = 100),rlnorm(1,0,10)))
  }
  if (!is.null(init) && length(init) != n.chains) stop("Initial values should be specified for every chain separately.")
  
  
  pe_bi[1,,] <- unlist(init) # intial values are the first row of the burn-in sample
  mu <- vector(length = p-1); t2 <- vector(length = p-1) #empty vectors for mu and tau squared
  var_params <- matrix(NA, n.iter, 2)
  
  if (length(which.MH != 0)) { # Only if MH is to be applied, the next errors hold
    if (length(which.MH) != length(MH.conditional) | length(which.MH) != length(MH.proposal) |
        length(which.MH) != length(MH.inital)) {
      stop("For every parameter that must be estimated by means of a Metropolis-Hastings step,
            a separate conditional distribution and a separate proposal distribution must be specified.")
    }
    if (!is.numeric(which.MH) | !is.numeric(MH.inital)) {
      stop("Input which.MH and MH.initial must be numeric vectors of length greater than zero, 
            containing the position of the parameter that has to be estimated by means of a 
            Metropolis-Hastings step.")
    }
    if (!is.list(MH.conditional) | !is.list(MH.proposal)) {
      stop("Both the conditional and the proposal must be specified as functions in a list, 
            even if there is only one function.")
    }
    # Set first "last-value" to initial value, and make object for acceptance rates
    MH.last <- matrix(MH.inital, n.chains, length(which.MH), byrow = TRUE)
    acctot <- matrix(0,n.chains,length(which.MH))
    
  }
  
  for (i in 2:burnin) { # for the number of burn-in samples specified
    
    if (i == 2) {cat(paste("\n", "Burn-in", "\n\n")); flush.console()}
    
    for (chain in 1:n.chains) { # for the number of chains
      
      b <- pe_bi[i-1,-p,chain] # specify initial regression coefficients and variance
      s2 <- pe_bi[i-1,p,chain] # based on the previous round of sampling
      
    #  if (ncol(x) == 2) b <- t(b) # if there are only two parameters, R sees a 
                                  # single value as a row instead of a column vector
      
      for (par in 1:(p-1)) { #for the number of parameters (-1 for the variance)
        
        if (par %in% which.MH) { #If M-H hastings is to be used, follow the following lines
          
          index <- which(which.MH %in% par) #Which parameter has to be used for MH (returns a single value)
    
          C <- MH.proposal[[index]]() #Draw a candidate value from the proposal density
          # Calculate the acceptance region by means of the user-specified conditional
          A <- MH.conditional[[index]](C, Y = y, X1 = x[,par], Xrest = x[,-par], B = b[-par], S2 = s2) - 
            MH.conditional[[index]](MH.last[chain, index], Y = y, X1 = x[,par], Xrest = x[,-par], B = b[-par], S2 = s2)
          
          # If the log of the random draw is smaller than the log of the acceptance region (which
          # works trough user-specified log-transformed conditonal input arguments) (and the minus
          # sign instead of a division in A), replace the last value of that chain by the candidate
          if (log(runif(1)) < A) {
            MH.last[chain, index] <- C 
          }
          b[par] <- MH.last[chain, index] #and update the parameter value by the MH.value
        }
        
        else {
          # update every variable on a variable by variable basis, as well as its variance
          mu[par]    <- ((sum(x[,par]*(y - x[,-par, drop=F]%*%b[-par]))/s2)+(m00[par]/t2_00[par]))/((sum(x[,par]^2)/s2) + (1/t2_00[par]))
          t2[par]    <- 1 / (sum(x[,par]^2)/s2 + 1/t2_00[par])
          # and calculate its regression coefficient
          b[par]     <- rnorm(1,mu[par],sqrt(t2[par]))
        }
      }
      
      # for every chain and every iteration, calculate the parameters of the gamma function
      A1 <- n/2 + A00; B1 <- sum((y - x%*%b)^2)/2 + B00
      # and sample the variance itself
      s2 <- 1 / rgamma(1,A1,B1)
      # and put the values in the output object for burn-in samples(might be useful to check convergence)
      pe_bi[i,-p,chain] <- b; pe_bi[i,p,chain] <- s2
    }
    setTxtProgressBar(pb_bi, value = i) #update the progress bar
  }
  
  for (i in 1:n.iter) { # for the number of iterations samples specified
    
    if (i == 1) {cat(paste("\n", "\n", "Sampling","\n", "\n")); flush.console()}
    
    for (chain in 1:n.chains) { # for the number of chains
      
  #    if (ncol(x) == 2) b <- t(b) # if there are only two parameters, R sees a 
      # single value as a row instead of a column vector
      
      for (par in 1:(p-1)) { #for the number of parameters (-1 for the variance)
        
        if (par %in% which.MH) { # If M-H hastings is to be used, for the specific paramter
                                 # follow the following lines
          index <- which(which.MH %in% par) #Which parameter has to be used for MH (returns a single value)
          
          C <- MH.proposal[[index]]() #Draw a candidate value from the proposal density
          # Calculate the acceptance region by means of the user-specified conditional
          A <- MH.conditional[[index]](C, Y = y, X1 = x[,par], Xrest = x[,-par], B = b[-par], S2 = s2) - 
            MH.conditional[[index]](MH.last[chain, index], Y = y, X1 = x[,par], Xrest = x[,-par], B = b[-par], S2 = s2)
          
          # If the log of the random draw is smaller than the log of the acceptance region (which
          # works trough user-specified log-transformed conditonal input arguments) (and the minus
          # sign instead of a division in A), replace the last value of that chain by the candidate
          # and increase the accepted total by 1
          if (log(runif(1)) < A) {
            MH.last[chain, index] <- C 
            acctot[chain, index] <- acctot[chain,index] + 1
          }
          b[par] <- MH.last[chain, index]
        }
        
        else {
          # update every variable on a variable by variable basis, as well as its variance
          mu[par]    <- ((sum(x[,par]*(y - x[,-par, drop=F]%*%b[-par]))/s2)+(m00[par]/t2_00[par]))/((sum(x[,par]^2)/s2) + (1/t2_00[par]))
          t2[par]    <- 1 / (sum(x[,par]^2)/s2 + 1/t2_00[par])
          # and calculate its regression coefficient
          b[par]     <- rnorm(1,mu[par],sqrt(t2[par]))
        }
      }
      
      # for every chain and every iteration, calculate the parameters of the gamma function
      A1 <- n/2 + A00; B1 <- sum((y - x%*%b)^2)/2 + B00; var_params[i,] <- c(A1,B1)
      # and sample the variance itself
      s2 <- 1 / rgamma(1,A1,B1)
      # and put the values in the output object for the sampled values
      pe_ni[i,-p,chain] <- b; pe_ni[i,p,chain] <- s2
    }
    setTxtProgressBar(pb_ni, value = i) #set progress bar for sampling
  }
  
  cat("\n"); flush.console()
  
  b_means <- rowMeans(colMeans(pe_ni))[-p]
  s_mean <- rowMeans(colMeans(pe_ni))[p]
  
  ## DIC = 2*Dbar - Dhat
  if (isTRUE(DIC.samples)) {
    D_hat <- -2 * (- (n/2) * log(2*pi) - n*log(sqrt(s_mean)) - (1/(2*s_mean)) * sum((y - x%*%b_means)^2))
    Dev <- vector("numeric", length = n.iter*n.chains)
    
    for (i in 1:n.iter) {
      for (k in 1:n.chains) {
        Dev[(i-1)*n.chains + k] <- -(n/2) * log(2*pi) - n*log(sqrt(pe_ni[i,p,k])) - (1/(2*pe_ni[i,p,k])) * sum((y - x%*%pe_ni[i,-p,k])^2)
      }
    }
    D_bar <- mean(-2*Dev)
    DIC <- D_hat + 2*(D_bar - D_hat)
  }
  
  autocor <- array(NA, dim = c(30,p,n.chains), dimnames = a_names) #create empty array for autocorrelations (30 max)
  for (i in 1:30) { #maximum of 30 autocorrelations
    for (par in 1:p) { #for every parameter
      for (k in 1:n.chains) { #and every chain
        autocor[i,par,k] <- cor(pe_ni[1:(n.iter-29),par,k], pe_ni[i:(n.iter-30+i),par,k])}}} #calculate autocor and end loops
  
  if (isTRUE(trace.plot)) {
    
    par.initial <- par()$mfrow #save initial graphical settings
    
    if (p <= 4)               par(mfrow = c(2,2)) #specify new settings depending on 
    else if (p > 4 && p <= 6) par(mfrow = c(3,2)) #the number of parameters
    else                      par(mfrow = c(3,3))
    
    for (par in 1:p) { #create trace plot for every parameter
      if (par == 1) {cat(paste("\n", "Convergence plots", "\n\n")); flush.console()} #print Convergence Plots to console
      plot(pe_ni[,par,1], type = "l", col = suppressWarnings(brewer.pal(n.chains, "Set1"))[1],
           main = colnames(pe_ni[,,1])[par], ylab="Estimate", xpd = FALSE)
      
      if (n.chains > 1) {       #if the number of chains is larger than 1
        for (k in 2:n.chains) { #add a line for every chain
          lines(pe_ni[,par,k], col = suppressWarnings(brewer.pal(n.chains, "Set1"))[k], xpd = FALSE)
        }
      }
      setTxtProgressBar(pb_cp, value = par) #and update the progress bar
    }
    par(mfrow = par.initial) #reset initial settings
    cat("\n"); flush.console()
  }
  if (isTRUE(autocor.plot)) {
    
    par.initial <- par()$mfrow #store initial graphical settings
    
    if (p <= 4)               par(mfrow = c(2,2)) # update graphical settings based on 
    else if (p > 4 && p <= 6) par(mfrow = c(3,2)) # the number of parameters
    else                      par(mfrow = c(3,3))
    
    for (par in 1:p) { #for every parameter
      if (par == 1) {cat(paste("\n", "Autocorrelation plots", "\n\n")); flush.console()} #print autocorrelation plots to console
      for (k in 1:n.chains) { #and every chain
        barplot(autocor[,par,k], space = 0, #create a barplot with the autocorrelations
                col = "forestgreen", main = paste0("Autocorrelation ", colnames(pe_ni)[par], ", chain ", k),
                ylab = "Autocorrelation", ylim = c(-1,1), axes = TRUE, axisnames = TRUE)
        axis(side = 1, at = seq(0,30,by=5), labels = seq(0,30,by=5)) #set an axis
        setTxtProgressBar(pb_ap, value = par) #and update the progressbar
      }
    }
    par(mfrow = par.initial) #and restore inital settings
  }
  
  param_est <- apply(pe_ni, 2, mean)
  param_sd  <- apply(pe_ni, 2, sd)
  param_ci  <- t(apply(pe_ni, 2, function(x) quantile(x, c(.025,.975))))
  param_naive <- apply(pe_ni, 2, function(x) {sd(x)/sqrt(n.chains*n.iter)})
  
  if (n.chains > 1) {
    s2j <- (1 / (n.iter-1)) * apply(pe_ni, 2, function(x) colSums((x - colMeans(x))^2))
    theta <- apply(pe_ni,2,mean)
    W <- (1/n.chains) * colSums(s2j)
    B <- (n.iter/(n.chains-1)) * (rowSums((apply(pe_ni,3,colMeans) - theta)^2))
  
    vartotal <- ((n.iter-1)/n.iter)*W + (1/n.iter) * B
    gelman_rubin <- sqrt(vartotal/W)
  }
  
  acc_rate <- rep(1,p)
  
  if (!is.null(which.MH)) {acc_rate[which.MH] <- colMeans(acctot)/n.iter}
  
  if (n.chains == 1) {
    summary_out <- round(data.frame(param_est,param_sd,param_ci,param_naive,acc_rate), digits = 3)
    names(summary_out) <- c("Mean", "SD", "2.5%", "97.5%", "Naive SE", "Acceptance rate")}
  else {
    summary_out <- round(data.frame(param_est,param_sd,param_ci,param_naive,gelman_rubin,acc_rate), digits = 3)
    names(summary_out) <- c("Mean", "SD", "2.5%", "97.5%", "Naive SE", "Gelman-Rubin", "Acceptance rate")}
  
  if (isTRUE(DIC.samples)) {
    output <- list(summary_out = summary_out, autocorrelation = autocor, posterior = pe_ni, DIC = DIC, D_bar = D_bar, 
                   n.iter = n.iter, n.burnin = burnin, n.chains = n.chains, call_object = call_object)}
  else {
    output <- list(summary_out = summary_out, autocorrelation = autocor, posterior = pe_ni, 
                 n.iter = n.iter, n.burnin = burnin, n.chains = n.chains, call_object = call_object)}
  
  class(output) <- "GibbsOut"
  
  return(output)  
}

############################################################################
## Third part, print function
############################################################################

print.GibbsOut <- function(x,...) {    # To get the output nice and organised
  mc <- x[["call_object"]]             # without printing a very large dataframe
  formula <- mc$form                   # with all sampled values, I decided to write
  summary_stats <- x[["summary_out"]]  # and additional little print function
  DIC <- x[["DIC"]]                    # so that a summary of the output if printed
  Dbar <- x[["D_bar"]]                 # to the console, but that all objects can be
  n.iter <- x[["n.iter"]]              # found in the output object itself
  burnin <- x[["n.burnin"]]
  chains <- x[["n.chains"]]
  cat(paste0("\n",
             "Bayesian linear regression", "\n\n",
             "Formula: ", paste(as.character(as.formula(formula))[c(2,1,3)], collapse = " "), "\n\n",
             "Iterations per chain = ", burnin + 1,":", burnin + n.iter, "\n",
             "Number of chains = ", chains, "\n",
             "Iterations per chains = ", n.iter, "\n",
             ifelse(!is.null(DIC), paste0(
             "DIC = ", round(DIC,3), "\n",
             "Dbar = ", round(Dbar,3), "\n\n"),""), 
             "Parameter estimates", "\n\n"))
  print(summary_stats)
}

############################################################################
## Fourth part, first / overall analyses
############################################################################

## Initial values must be specified in a list of lists, containing one list per chain
init1 <- list(b0 = 0, b1 = 0, b2 = 0, b3 = 0, sigma_squared = 1)
init2 <- list(b0 = 100, b1 = 100, b2 = 100, b3 = 100, sigma_squared = 100)
inits <- list(init1, init2)

## When an informative prior is used, it must be specied in a list of vectors. The first vector
## must contain the means the regression coefficients (including intercepts),
## the second vector should contain the corresponding standard deviations. The third 
## vector should contain the prior alpha and beta for the variance parameter.

means <- c(0, 6, 1, 1)
vars  <- c(2, .5, .1, 1)
sigma_squared <- c(100, 20)

prior <- list(means, vars, sigma_squared)

## When a metropolis hastings step is to be used, it has to be of the form
## function(P(roposal value), Y, X, B, S2)

proposal_density <- list(function() rnorm(1, 0.028, .13))

## Make sure to log transform the conditional density, and make sure that the 
## conditional has input arguments for the proposal value, the DV, the variable of 
## interest in this step (X1), all other IV's (Xrest) the regression coefficients 
## of the other variables (B), and the variance (S2). All other arguments 
## must be specified in the function itself
conditional_density <- list(function(P, Y, X1, Xrest, B, S2, m0 = 0, s02 = 10000, v0 = 2) 
  {(-P^2*(sum(X1^2))/(2*S2) + P*(sum(X1*(Y - Xrest%*%B)))/S2) +
    log((1+((P-m0)^2)/(v0*s02))^((-(v0+1)/2)))})


## Note that in "which.MH", the value 1 is indicative for the intercept, and the value
## "number of regression coefficients + 2" would be indicative for the residual variance
out1 <- bayes_lm(form = watts ~ nact + speed + distance, data = dat, prior = "uninf",
                which.MH = 4, MH.inital = 1, MH.conditional = conditional_density,
                MH.proposal = proposal_density, 
                DIC.samples = TRUE,
                burnin = 1000, n.iter = 10000, n.chains = 3, 
                seed = 12345, trace.plot = TRUE, autocor.plot = TRUE)

out2 <- bayes_lm(form = watts ~ nact + speed, data = dat, prior = "uninf",
                 DIC.samples = TRUE,
                 burnin = 1000, n.iter = 10000, n.chains = 3, 
                 seed = 12345, trace.plot = TRUE, autocor.plot = TRUE)

out3 <- bayes_lm(form = watts ~ speed, data = dat, prior = "uninf",
                 DIC.samples = TRUE,
                 burnin = 1000, n.iter = 10000, n.chains = 3, 
                 seed = 12345, trace.plot = TRUE, autocor.plot = TRUE)

############################################################################
## Fifth part, Bayes Factor function
############################################################################

BF <- function(GibbsOut, hypos = NULL, n.constraints = NULL, region_of_equivalence = 0.5, n.samples = 100000) {
  
  # This Bayes Factor function only works with the output of my Gibbs sampling function
  if(class(GibbsOut) != "GibbsOut") stop("Input must be of class GibbsOut.")
  
  x <- GibbsOut          # get the gibbs output from the sampling function
  mc <- x$call_object    # extract the function call from the gibbs output
  d <- eval(mc$data)     # now we have the data from the gibbs output as well
  f <- mc$form           # and we additionally have the formula
  ud <- model.frame(f,d) # so that we can obtain the actual data.frame from the gibbs function
  
  posterior <- x$posterior    # posterior samples so that we can calculate the means
  n.params <- ncol(posterior) # extract the number of parameters
  array_to_matrix <- apply(posterior,2,rbind) # and make one long matrix from the initial array containg all chains
  reg_est <- array_to_matrix[,-n.params] # additionally, remove the variance, because this will not be tested
  
  # specify which parameters have to be included (those that are mentioned in the hypotheses)
  test_which <- colnames(reg_est) %in% unlist(str_split(hypos, pattern = "([^[:alnum:]])"))
  # and only keep those estimates from the posterior sampled values
  reg_est <- reg_est[,test_which]
  
  # make sure that the data is standardized, if not, stop the function
  if(!isTRUE(all.equal(sapply(ud[,-1],mean), rep(0,n.params-2),check.attributes=F))) {
    stop("Data must be standardized in order to meaningfully compare regression coefficients.")}
  if(!isTRUE(all.equal(sapply(ud[,-1],var), rep(1,n.params-2),check.attributes=F))) {
    stop("Data must be standardized in order to meaningfully compare regression coefficients.")}
  
  varcov <- var(reg_est)     # posterior covariance between the parameter estimates
  means <- colMeans(reg_est) # posterior means of the parameter estimates
  
  # Now, sample values for the fit and complexity of the hypothesis
  # this is not yet dependent on the hypotheses itself
  fit_samples <- data.frame(MASS::mvrnorm(n.samples,means,varcov)) 
  complexity_samples <- data.frame(MASS::mvrnorm(n.samples,rep(0,length(means)),varcov))
  
  # slightly adjust the hypothesis if it does not contain spaces between the arguments
  create_spaces <- gsub("([^[:alnum:][:blank:]])"," \\1 ",hypos)
  # separate the hypotheses, so that every hypothesis is a single element in a vector of hypotheses
  sep_hypos <- strsplit(create_spaces, split = " ; ")
  # remove the "\n" statement from the hypothesis, which occurs due to specifying each hypothesis on a new line
  hypo <- str_trim(sapply(sep_hypos, function(x) gsub("\n","",x)))
  # create an empty matrix for the samples that are in accordance with the hypothesis (fit) and the prior (complexity)
  fit <- complexity <- matrix(NA,n.samples,length(hypo))
  # get the standard deviations from the covariance matrix
  sds <- sqrt(diag(varcov))
  # Now first, attach the samples fit, that is, the sampled values based on the estimated means
  attach(fit_samples, warn.conflicts = F)     # Attach is needed so that the eval(parse(text = ...)) command works
  for (i in 1:length(hypo)) {                 # Then, for all hypotheses, evaluate the following statements
    if(grepl("=", hypo[i])) {                 # If the hypothesis contains an equality constraint
       vars <- unlist(strsplit(hypo[i], " ")) # extract the variables from this specific hypothesis
       pool_which <- names(sds) %in% vars     # and check which sd's belong to these elements
       pooled_sd <- mean(sds[pool_which])     # then, pool the sds to possibly obtain a region of equivalence in terms of the pooled sd
       hypo[i] <- gsub("=", "-", hypo[i])     # Then, replace the equality by a minus sign
                                              # And in the next step, add that it concerns the absolute difference between the hypos
       hypo[i] <- paste("abs(",hypo[i],") < ", region_of_equivalence, "*" , pooled_sd) # that must be smaller than the region of equivalence
    }                                         # To wrap up, we test whether the absolute difference between the sampled values
    fit[,i] <- eval(parse(text = hypo[i]))    # is smaller than the region of equivalence, which is evaluated on this line
  }
  detach(fit_samples)
  for (i in 1:length(hypo)) {                    # Once again, loop through every hypothesis
    if(is.null(n.constraints)) j <- str_count(hypo[i], pattern = "\\<|\\>|\\=")
    else j <- n.constraints[i]                   # number of independent constraints
    frac <- j/nrow(ud)                           # fraction is j divided by sample size
    if(grepl("=", hypo[i])) {                    # Check whether there is an equality constrain
      vars <- unlist(strsplit(hypo[i], " "))     # Split the hypothesis at the spaces
      pool_which <- names(sds) %in% vars         # Check which variables the current hypothesis concerns
      pooled_sd <- mean(sds[pool_which])         # And pool the standard deviations of the regression coefficients
      
      hypo[i] <- gsub("=", "-", hypo[i])         # Now, replace the equality by a minus sign
                                                 # And add that we are interested in the absolute difference
      hypo[i] <- paste("abs(",hypo[i],") < ", region_of_equivalent, "*", pooled_sd)
    }
    complexity_samples <- data.frame(MASS::mvrnorm(n.samples, rep(0, length(means)), varcov / frac))
    attach(complexity_samples, warn.conflicts = F) # attach complexity to make sure that the hypotheses are calculated on this object
    complexity[,i] <- eval(parse(text = hypo[i])) # And compute the samples that are in accordance with the hypothesis
  } 
  detach(complexity_samples)               # We do not need to look in the complexity samples anymore, so detach
  
  output <- data.frame(Fit = colMeans(fit), Com = colMeans(complexity), BFu = colMeans(fit)/colMeans(complexity),
                       BFc = (colMeans(fit)/colMeans(complexity)) / ((1-colMeans(fit))/(1-colMeans(complexity))),
                       PMP = (colMeans(fit) / colMeans(complexity)) / sum(colMeans(fit)/colMeans(complexity)))
  rownames(output) <- paste0("H",1:length(hypo))
  return(list(BayesFactor = output))
}

############################################################################
## Sixth part, calculate Bayes Factors
############################################################################

scaled_data <- as.data.frame(lapply(model.frame(watts ~ nact + speed + distance, data = dat),scale))

BF_out <- bayes_lm(form = watts ~ nact + speed + distance, data = scaled_data, prior = "uninf",
                         burnin = 1000, n.iter = 10000, n.chains = 6, 
                         seed = 12345, trace.plot = T, autocor.plot = T, DIC.samples = T)

hypos1 <- "speed > 0 & nact > 0;
           nact > speed"

# The region of equivalence is stated in terms of the (pooled) standard error(s) of the
# parameters of interest. 
BFs1 <- BF(BF_out, hypos1, region_of_equivalence = 2, n.samples = 1000000)

hypos2 <- "speed > 0 & nact > 0;
           nact > speed;
           speed & nact & distance"

BFs2 <- BF(BF_out, hypos2, n.constraints = c(2,1,1), region_of_equivalence = 2, n.samples = 1000000)

############################################################################
## Seventh part, calculate posterior predictive p-value
############################################################################

# First, extract the data used in the analysis
used_data <- model.matrix(formula(out1$call_object$form), eval(out1$call_object$data))
# Then, extract the used outcome variable
obs_y <- model.frame(formula(out1$call_object$form), eval(out1$call_object$data))[,1]
# Set a seed for reproducibility
set.seed(100)
## For the simulated residuals: Residuals = observed score - predicted score
## Observed score = predicted score + random draw from a N(0,sigma^2)
## distribution, so, the residuals are equal to the random draw from the N(0,sigma^2)
sim_res <- apply(matrix(out1$posterior[,5,],ncol=1), 1, function(x) rnorm(obs_y,0,sqrt(x)))
# Calculate the observed residuals for all iterations of sampled values
obs_res <- apply(apply(out1$posterior,2,rbind),1,function(x) obs_y - used_data %*% x[-5])
# Calculate proportion of the skewness of the  residuals of simulated datasets that
# is larger than the skewness of the observed residuals over the sampled parameter
# estimates.
skewness_p <- mean(abs(apply(sim_res,2,skewness)) > abs(apply(obs_res,2,skewness)))

cor_p <- mean(apply(sim_res,2, function(x) cor(qnorm(ppoints(x)),sort(x))) >
     apply(obs_res,2, function(x) cor(qnorm(ppoints(x)),sort(x))))



############################################################################
## Eighth part, Bayesian Sequential Updating
############################################################################
## In the meantime of performing these analyses, I have made some additional
## bikerides, that were not included in the previous analysis yet. Since
## from a Bayesian viewpoint, incorporating previous knowledge is easy,
## the posterior from the previous analysis will be used in the current
## analysis. 

up_dat <- dat_uncentered
up_dat[,c(2,3,5,6)] <- as.data.frame(lapply(up_dat[,c(2,3,5,6)], FUN = function(x) x-mean(x, na.rm = TRUE)))


## linear model with the data, regress watts on number of activities in the past four 
## weeks, the distance of the ride, and the average speed of the ride. Only the latter
## subset of the data
up_out <- bayes_lm(form = watts ~ nact + speed + distance, data = up_dat[89:107,], prior = "uninf",
                   DIC.samples = TRUE,
                   burnin = 1000, n.iter = 10000, n.chains = 3, 
                   seed = 12345, trace.plot = TRUE, autocor.plot = TRUE)

## linear model with the complete data, regress watts on number of activities in the past four 
## weeks, the distance of the ride, and the average speed of the ride.
full_out <- bayes_lm(form = watts ~ nact + speed + distance, data = up_dat, prior = "uninf",
                   DIC.samples = TRUE,
                   burnin = 1000, n.iter = 10000, n.chains = 3, 
                   seed = 12345, trace.plot = TRUE, autocor.plot = TRUE)


