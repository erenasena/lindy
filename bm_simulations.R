### Survival and hazard functions of hitting times 

# The necessary packages
#install.packages("survminer")
#install.packages("survival")
#install.packages("dplyr")

# The necessary libraries
library(parallel)
library(survminer)
library(survival)
library(dplyr)
library(dynpred)

### The functions that prepare the data for the hypothesis tests, depending on what is tested

## Computes the differences between successive hazard rates to get the mean difference
diff <- function(x){ # the input is a vector of hazard rates 
  diff <- numeric(length = length(x) - 1)
  for(i in 2:length(x)){
    diff[i-1] <- x[i] - x[i-1]
  }
  return(diff)
}

## Tracks the changes in the hazard rates in terms of increases and decreases 
## Can be used both to get the total number of increases / decreases or successive increases / decreases
change <- function(x){ # the input is a vector of hazard rates
                       # returns a vector of 1s and -1s 
                  
  change <- numeric(length(x) - 1)
  for(i in 2:length(x)){
    if(x[i] - x[i-1] > 0){
      change[i-1] <- 1 # if the hazard rate increased, we get a 1  
    } else if(x[i] - x[i-1] < 0){
      change[i-1] <- -1 # if it decreased, we get a -1
    } else if(x[i] - x[i-1] == 0){
      change[i-1] <- 0 # if it remained the same, which it never does, we get a 0 
    }
  }
  return(change)
}

## If the goal is to get the successive increases / decreases, this function will do it
reps <- function(x){ # the input is the output vector of the change function  
  rep_vals <- rle(x)$values # just the values (i.e., 1 and -1)
  rep_lengths <- rle(x)$lengths # how many times they are repeated
  data <- data.frame('Values' = rep_vals, 'Frequency' = rep_lengths) # all in a nice form 
  plus <- which(data$Values == 1) # indexing the values 
  minus <- which(data$Values == -1)
  inc <- sum(data$Frequency[plus]) # the total number of successive increases 
  dec <- sum(data$Frequency[minus]) # the total number of successive decreases 
  max_inc <- max(data$Frequency[plus]) # the maximum number of successive increases
  max_dec <- max(data$Frequency[minus]) # the maximum number of successive decreases 
  list = list('Data' = data, 'Decreases' = max_dec, 'Increases' = max_inc)
  return(list$Decreases)
}

### Generate trial data from the exponential, Pareto and Weibull distributions 
set.seed(1)
pareto <- VGAM::rpareto(n = 10000, scale = 1, shape = 1.5)

set.seed(1)
exps <- rexp(n = 10000, rate = 0.0001)

set.seed(1)
weibull <- rweibull(n = 10000, shape = 1.1)

## Explicit hazard functions 
time <- c(1:1000)

# The Weibull distribution 
wh <- function(a, b, t){
  h <- (b/a)*((t/a)^(b-1))
  return(h)
}

w <- wh(a = 1, b = 4, t = time)
plot(w, type = 'l', xlim = c(min(time), max(time)), ylim = c(min(w), max(w)))

# The exponential distirbution
eh <- function(x, b){
  e <- numeric(length = length(x))
  for(i in x){
    x[i] <- 1/b
  }
  return(e)
}

e <- eh(x = time, b = 1)
plot(e, type = 'l')

# Pareto hazard manually 
par <- function(x, a, b){
  fx <- VGAM::dpareto(x, shape = a, scale = b)
  dx <- VGAM::ppareto(x, shape = a, scale = b)
  sx <- 1-dx
  hx <- fx/sx
  return(hx)
}

p <- par(x = time, a = 2, b = 1)
plot(p, type = 'l')

# Pareto hazard function from Eliazar (2017)
ph <- function(p, t){
  h <- ((1+p) / p)*(1/t)
  return(h)
}

p <- ph(p = 1, t = time)
plot(p, type = 'l', main = 'eli')

## Calculating the hazards with the package 
survival_object <- Surv(time = pareto) 
survival_fit <- survfit(survival_object ~ 1)
hazard_fit <- bshazard::bshazard(survival_object ~ 1)
hazard_time <- hazard_fit$time
hazard_rates <- hazard_fit$hazard
new_plot <- plot(x = hazard_time, y = hazard_rates, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 100), ylim = c(min(hazard_rates), max(hazard_rates)))

## Hypothesis test on the mean difference 

# Resampling the hazard rates, getting their means, getting their differences, finding the mean 
# difference, and replicating this process 10,000 times to get a distribution of mean differences
mean_diff <- mean(diff(x = hazard_rates))
diff_dist <- replicate(10000, mean(diff(x = sample(hazard_rates, length(hazard_rates), FALSE))))
hist(diff_dist, breaks = 100, xlab = 'Mean change in successive hazard rates', main = 'Distribution of mean differences between sucessive hazard rates') # a histogram showing the distribution 
abline(v = mean_diff, col = 'red', lwd = 2) # drawing a line to locate our observed mean 

# One sided, the alternative is <
sum(diff_dist < mean_diff) / 10000 

# One sided, the alternative is > 
sum(diff_dist > mean_diff) / 10000

# Two tailed tests 
(sum(diff_dist < mean_diff) + sum(diff_dist > abs(mean_diff))) / 10000 # two tailed for smaller 
sum(abs(diff_dist) > mean_diff) / 10000 # two tailed for larger 

## Hypothesis test on the total number of decreases 

# Resampling for the sums 
sum_minus <- sum(change(hazard_rates) == -1) # the observed total number of decreases
sum_dist <- replicate(10000, sum(change(x = sample(hazard_rates, length(hazard_rates), FALSE)) == -1)) 
hist(sum_dist, breaks = 100, xlab = 'Total number of decreases')  
abline(v = sum_minus, col = 'red', lwd = 2)

# One sided, the proportion of total decreases larger than the observed sum 
sum(sum_dist > sum_minus) / 10000

## Hypothesis test on the difference between the total number of increases and decreases 

# Resampling 
RNGkind("L'Ecuyer-CMRG") # this is the random number generator needed in parallel processing 
detectCores() # tells you the number of cores your computer can use for the simulations 

sum_diff <- sum(change(hazard_rates) == 1) - sum(change(hazard_rates) == -1) # the observed difference

sum_diff_fun <- function(x){
  resampled <- sample(hazard_rates, length(hazard_rates), FALSE)
  changes <- change(x = resampled)
  plus <- sum(changes == 1)
  minus <- sum(changes == -1)
  return(plus - minus)
}

sum_diff_dist <- mclapply(1:10000, sum_diff_fun, mc.cores = 8, mc.set.seed = TRUE)
sum_diff_dist <- unlist(sum_diff_dist)
hist(sum_diff_dist, breaks = 500, xlab = 'Increases - decreases')  
abline(v = sum_diff, col = 'red', lwd = 2)

## Hypothesis test on the maximum number of successive increases or decreases 
observed_max <- reps(change(hazard_rates))
max_sucs_dist <- replicate(10000, reps(change(x = sample(hazard_rates, length(hazard_rates)))))
hist(max_sucs_dist, breaks = 100, xlab = 'Max number of successive decreases')  
abline(v = observed_max, col = 'red', lwd = 2)







### The BM functions 

## Arithmetic Brownian Motion (resamples until the value is below the boundary)  
my_abm <- function(nsim, t0, t, n, X0, mu, sigma, L, R){ 
  
  dt <- t/n 
  time <- seq(from = t0, to = t, by = dt)   
  initial <- X0  
  X <- matrix(nrow = nsim, ncol = length(time)) 
  X[1:nsim, 1] <- X0
  event_time <- numeric(length = nsim) 
  event_status <- numeric(length = nsim) 
  
  for(i in 1:nrow(X)){
    for(j in 2:length(time)){
      X[i,j] <- X0 + mu * dt + sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1)
      
      while(X[i,j] > R){ # simulate data as long as it's above the reflecting barrier 
        X[i,j] <- X0 + mu * dt + sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1)
        
        if(X[i,j] <= R) { # if we are below or equal, the value stays the same
        X[i,j] <- X[i,j]
        }
      }
      
      if(X[i,j] > L & j < ncol(X)){ 
        X0 <- X[i,j] 
        
        } else if(X[i,j] > L & j == ncol(X)){
        X0 <- initial 
        event_time[i] <- ncol(X) 
        event_status[i] <- 0 
        
        } else if(X[i,j] <= L){ 
        X[i,j] <- L
        X0 <- initial 
        event_time[i] <- j 
        event_status[i] <- 1
        break
      }
    }
    values <- list("Values" = X, "Event time" = event_time, "Event status" = event_status)
  }
    return(values)
}


## Geometric Brownian Motion (the new value is the barrier value until it is below)
my_gbm <- function(nsim, t0, t, n, X0, mu, sigma, L, R){ 
  dt <- t/n 
  sig2 <- sigma^2 
  time <- seq(from = t0, to = t, by = dt)  
  
  initial <- X0  
  X <- matrix(nrow = nsim, ncol = length(time)) 
  X[1:nsim, 1] <- X0 
  event_time <- numeric(length = nsim)
  event_status <- numeric(length = nsim)
  
  for(i in 1:nrow(X)){
    for(j in 2:length(time)){
      X[i,j] <- X0 * exp(((mu - 0.5 * sig2) * dt) + (sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1))) # the formula for the next value of BM
      
      if(X[i,j] > R){ # if above the barrier, the new value is the barrier level 
        X[i,j] <- R

        } else if(X[i,j] <= R) { # otherwise it's the simulated value 
          X[i,j] <- X[i,j]
          }
      
      if(X[i,j] > L & j < ncol(X)){ 
        X0 <- X[i,j] 
        
      } else if(X[i,j] > L & j == ncol(X)){
        X0 <- initial 
        event_time[i] <- ncol(X) 
        event_status[i] <- 0 
        
      } else if(X[i,j] <= L){ 
        X0 <- initial 
        event_time[i] <- j 
        event_status[i] <- 1 
        break
        
      }
    }
    values <- list("Values" = X, "Event time" = event_time, "Event status" = event_status)
  }
  return(values)
}

### Simulations
RNGkind("L'Ecuyer-CMRG") # this is the random number generator needed in parallel processing 
detectCores() # tells you the number of cores your computer can use for the simulations 

## Storing the values: These two functions will allow us to nealty store the resutls. It is easier to 
## define them before running the simulations. I wrote them because the results of parallel simulations are a mess. 

# The BM values 
values <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) == "Event time")]
  u <- u[-which(names(u) == "Event status")]  
  m <- matrix(data = u, nrow = nsim, ncol = length(0:n), byrow = TRUE)
  df <- data.frame(x = t(m), row.names = paste0("Time", 0:n, ""))
  colnames(x = df) <- paste0("Sim", 1:nsim, "")
  store <- list(m, df)
  return(store)
}

# The hitting times 
times <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "Event time")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("Sim", 1:nsim, ""))
  colnames(x = df) <- "Hitting time"
  store <- list(m, df)
  return(store)
}

# The event / censoring times 
events <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "Event status")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("Sim", 1:nsim, ""))
  colnames(x = df) <- "Event"
  store <- list(m, df)
  return(store)
}

## Parallel runs 
f <- function(i){ # specify the desired function and parameter values here
  my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90, R = 10000000) 
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 1, mu = -1, sigma = 1, L = 0.65, R = 1.01)
}

set.seed(1)
res <- mclapply(X = 1:10000, f, mc.cores = 8, mc.set.seed = TRUE)

v <- values(x = res, nsim = 10000, n = 1000) # indexing the BM values 
m_val <- v[[1]] # BM values in a matrix (goes into the plotting function)
df_val <- v[[2]] # BM values in a data frame

t <- times(x = res, nsim = 10000, n = 1000) # indexing the hitting times 
m_times <- t[[1]] # in a matrix (for histograms)
df_times <- t[[2]] # in a data frame 

e <- events(x = res, nsim = 10000, n = 1000)
m_event <- e[[1]] # in a matrix
df_event <- e[[2]]

# Histogram of hitting times
hist(m_times, breaks = 100, xlim = c(0, 1010), main = 'GBM with an absorbing barrier')
legend(x = "center", legend = c('mu = -1', 'sigma = 1', 'L = 90', 'R = -100'))

### Survival curves and fitting the model 
surv_data <- data.frame(Time = m_times, Event = m_event, row.names = paste0("Sim", 1:nrow(m_times), ""))
surv_object <- Surv(time = m_times, event = m_event) 
surv_fit <- survfit(surv_object ~ 1)
#surv_plot <- plot(x = surv_fit$time, y = surv_fit$surv, type = 'l', xlab = "Time")
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)
hazard <- haz_fit$hazard
hazard_plot <- plot(haz_fit$time, hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 1000), ylim = c(min(haz_fit$haz), max(haz_fit$haz)))