### Survival and hazard functions of hitting times 

# The necessary packages
#install.packages("survminer")
#install.packages("survival")
#install.packages("dplyr")

# The necessary libraries 
library(ggplot2)
library(reshape2)
library(parallel)
library(survminer)
library(survival)
library(dplyr)
library(dynpred)

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

## Geometric Brownian Motion (reflects back; the new value is the previous value)
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
      
      if(X[i,j] > R){ # if above the barrier, the value does not change 
        X[i,j] <- X[i,j-1]

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



### Visualizations
bmplot <- function(x, nsim, n, L, R, ylim, title){ # x is the matrix output of the BM functions, n is the number of simulations, t is the vector of time points, as in the BM functions 
  rownames(x) <- paste("sim", seq(nsim), sep = "") # the number of simulations / rows
  colnames(x) <- paste("time", seq(0:n), sep = "") # the number of time points - 0:100 at the moment / columns 
  dat <- as.data.frame(x) # creating the data frame for ggplot 
  dat$sim <- rownames(dat)
  mdat <- melt(dat, id.vars = "sim")
  mdat$time <- as.numeric(gsub("time", "", mdat$variable))
  
  p <- ggplot(data = mdat, mapping = aes(x = time, y = value, group = sim)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(angle = 0, size = 11, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
          axis.title.x = element_text(margin = margin(t = 10, b = 10))) +
    geom_line(size = 0.3, alpha = 1, aes(color = sim), show.legend = FALSE) + 
    ggtitle(title) + xlab("Time") + ylab("Value") + ylim(ylim) + 
    geom_hline(yintercept = L, color = "orange", size = 0.5) +
    geom_hline(yintercept = R, color = "green", size = 0.5)
  return(p)
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

#p <- bmplot(x = m_val, nsim = 10000, n = 100, L = 90, R = 75, ylim = c(min(m_val), max(m_val)), # Define the range of the y-axis  
            #title = "Brownian motion with an absorbing barrier")
#print(p)

# Histogram of hitting times
hist(m_times, breaks = 100, xlim = c(0, 1010), main = 'GBM with an absorbing barrier')
legend(x = "center", legend = c('mu = -1', 'sigma = 1', 'L = 90', 'R = -100'))

### Survival curve and fitting the model 
surv_data <- data.frame(Time = m_times, Event = m_event, row.names = paste0("Sim", 1:nrow(m_times), ""))
surv_object <- Surv(time = m_times, event = m_event) 
surv_fit <- survfit(surv_object ~ 1)
#surv_plot <- plot(x = surv_fit$time, y = surv_fit$surv, type = 'l', xlab = "Time")
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)
hazard <- haz_fit$hazard
hazard_plot <- plot(haz_fit$time, hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 1000), ylim = c(min(haz_fit$haz), max(haz_fit$haz)))

## A function that computes the difference between successive rates 
diff <- function(x){
  diff <- numeric(length = length(x) - 1)
  for(i in 2:length(x)){
    diff[i-1] <- x[i] - x[i-1]
  }
  return(diff)
}

diff_haz <- diff(x = hazard) # a vector of the differences between successive hazard rates 
observed <- mean(diff_haz)
median(diff_haz)
plot(diff_haz, type = 'l')

### Generate perfect data from the exponential, Pareto and Weibull distributions 
set.seed(1)
pareto <- VGAM::rpareto(n = 10000, scale = 1, shape = 1.5)

set.seed(1)
exps <- rexp(n = 10000, rate = 0.0001)

set.seed(1)
weibull <- rweibull(n = 10000, shape = 0.9)
exp_weibull <- function(t, lambda = 0.5, gamma = 1.5) lambda * gamma * t^(gamma - 1)

## Calculating the hazards 
new_surv_object <- Surv(time = pareto) 
new_surv_fit <- survfit(new_surv_object ~ 1)
new_fit <- bshazard::bshazard(new_surv_object ~ 1)
new_hazard <- new_fit$hazard
new_plot <- plot(new_fit$time, new_fit$hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 100), ylim = c(min(new_fit$haz), max(new_fit$haz)))
new_diff <- diff(x = new_hazard)
trial_observed <- mean(new_diff)

## Resampling the mean difference 
dist <- replicate(10000, mean(diff(x = sample(new_hazard, length(new_hazard), FALSE))))
hist(dist, breaks = 100)  
abline(v = trial_observed, col = 'red', lwd = 2) 

## Hypothesis tests

# One sided, the alternative is <
sum(dist < trial_observed) / 10000 # one sided, alternative is smaller 

# One sided, the alternative is > 
sum(dist > trial_observed) / 10000 # one sided, alternative is larger 

# Two tailed tests in different ways 
(sum(dist < trial_observed) + sum(dist > abs(trial_observed))) / 10000 # two tailed for smaller 
sum(abs(dist) > trial_observed) / 10000 # two tailed for larger 

## A function to get the number of successive decreases
change <- function(x){ # if the hazard rate increased, we get a 1, if it decreased, we get a -1 
  change <- numeric(length(x) - 1)
  for(i in 2:length(x)){
    if(x[i] - x[i-1] > 0){
      change[i-1] <- 1
      } else if(x[i] - x[i-1] < 0){
        change[i-1] <- -1
      } else if(x[i] - x[i-1] == 0){
        change[i-1] <- 0
      }
  }
  return(change)
}

changes <- change(x = new_hazard) # a vector of 1s and -1s 

reps <- function(x){
  rep_vals <- rle(x)$values # just the values (i.e., 1 and -1)
  rep_lengths <- rle(x)$lengths # how many times they are repeated
  data <- data.frame('Values' = rep_vals, 'Frequency' = rep_lengths) # all in a nice form 
  plus <- which(data$Values == 1) # indexing the values 
  minus <- which(data$Values == -1)
  inc <- sum(data$Frequency[plus]) # the total number of increases in succession 
  dec <- sum(data$Frequency[minus]) # the total number of decreases in succession 
  max_inc <- max(data$Frequency[plus]) # the number of first sucessions of 1 
  max_dec <- max(data$Frequency[minus]) # the number of first sucessions of -1 
  list = list('Data' = data, 'Decreases' = max_dec, 'Increases' = max_inc)
  return(list)
}

obs_rep <- reps(x = changes)

## Resampling for the repetitions 
rep_dist <- replicate(10000, reps(change(x = sample(new_hazard, length(new_hazard), FALSE))))
hist(rep_dist, breaks = 100)  
abline(v = obs_rep, col = 'red', lwd = 2)
sum(dist > obs_rep) / 10000

newhazards <- sample(new_hazard, length(new_hazard), FALSE)
plot(new_fit$time, newhazards, type = 'l')
newchanges <- change(x = newhazards)
newreps <- reps(x = newchanges)

