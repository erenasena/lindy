# The necessary libraries
library(dplyr) # don't know why 
library(parallel) # to run the parallel simulations 
library(ggplot2) # to plot the Brownian motion simulations
library(survival) # survival analysis
library(survminer) # to plot the survival curves 
library(dynpred) # to compute conditional survival 
library(qualityTools) # for the Q-Q plots 

### The Brownian motion functions 

## Arithmetic Brownian Motion (resamples until the value is below the reflecting barrier)  
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


## Geometric Brownian Motion (the new value equals the reflecting barrier until it is below it)
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

# A function to store the BM values in a nice form 
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

# A function to store the hitting times 
times <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "Event time")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("Sim", 1:nsim, ""))
  colnames(x = df) <- "Hitting time"
  store <- list(m, df)
  return(store)
}

# A function to store the event / censoring times 
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
  my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90, R = 11000000000) 
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 0, mu = 0, sigma = 1, L = -0.1, R = 10000000)
}

set.seed(1)
res <- mclapply(X = 1:8973, f, mc.cores = 8, mc.set.seed = TRUE)

v <- values(x = res, nsim = 8973, n = 1000) # indexing the BM values 
m_val <- v[[1]] # BM values in a matrix (goes into the plotting function)
df_val <- v[[2]] # BM values in a data frame

t <- times(x = res, nsim = 8973, n = 1000) # indexing the hitting times 
m_times <- t[[1]] # in a matrix (for histograms)
df_times <- t[[2]] # in a data frame 

e <- events(x = res, nsim = 8973, n = 1000)
m_event <- e[[1]] 
df_event <- e[[2]]

data <- cbind(df_times, df_event)

## Guessing the distribution 

# Histogram of the hitting times
hist(data$`Hitting time`, breaks = 10, xlim = c(0, 1010), main = 'GBM with an absorbing barrier')
legend(x = "center", legend = c('mu = -1', 'sigma = 1', 'L = 90', 'R = -100'))

# Descriptive statistics of the hitting times
quantile(data$`Hitting time`) # the sample quantiles
ext <- m_times[which(data$`Hitting time` > 250)] # the extreme values 
length(ext) / length(data$`Hitting time`) # what proportion of data are larger than a certain value 
sum(ext) / sum(data$`Hitting time`) # what proportion of the sum they make

# Q-Q Plot to check exponentiality (if concave, there might be heavy tailedness)
hittings <- sort(data$`Hitting time`) # sort the data
p <- ppoints(hittings, length(hittings)) # get the probabilities of the data
s <- quantile(x = hittings, p = p) # sample quantiles
q <- qexp(p = p) # exponential quantiles 
qqplot(x = s, y = q, xlab = "Sample quantiles", ylab = "Theoretical quantiles", 
       main = "Exponential QQ Plot")


qqPlot(x = data$`Hitting time`, y = "exponential", xlab = "Sample quantiles", 
       ylab = "Theoretical quantiles", main = "Exponential Q-Q Plot") # shorter way with confidence bands


# Zipf plot to check for power law decay 
hittings <- sort(data$`Hitting time`) # sort the data
fit <- ecdf(hittings) # empirical cumulative distribution function
s <- 1 - fit(hittings) # empirical survival function  
logs <- log(s) # the log of the survival function 
logx <- log(hittings) # the log of the sorted failure times 

plot(x = logx, y = logs, # if linear, indication of power law 
     xlab = "log(failure times)", ylab = "log(survival function)", 
     main = "The Zipf Plot") 

### Survival curves and fitting the model 
surv_data <- data.frame(Time = m_times, Event = m_event, row.names = paste0("Sim", 1:nrow(m_times), ""))
surv_object <- Surv(time = m_times, event = m_event) 
surv_fit <- survfit(surv_object ~ 1)
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)

survival <- surv_fit$surv
surv_time <- surv_fit$time

hazard <- haz_fit$hazard
haz_time <- haz_fit$time

hazard_plot <- plot(x = haz_time, y = hazard, xlab = 'Time', ylab = 'Hazard Rate', type = 'l', 
                    xlim = c(min(haz_time), max(haz_time)), ylim = c(min(haz_fit$haz), max(haz_fit$haz)))

## Conditional survival
fit <- Fwindow(object = surv_fit, width = 10, variance = TRUE, conf.level = 0.95)
con_time <- fit$time # the calculated times
con_death <- fit$Fw # conditional death 
con_surv <- 1 - con_death # conditional survival; they are mirror images
plot(x = con_time, y = con_surv, type = 'l', col = 'green', xlab = 'Time', 
     ylab = 'Probability', main = 'Conditional survival and death over time',
     ylim = c(0, 1), xlim = c(min(con_time), max(con_time)))
lines(con_time, con_death, col = 'red') # change the limit of the y-axis to c(0, 1) to see this 
cond <- data.frame(con_time, con_surv, con_death)
colnames(cond) <- c('Time', 'Conditional Survival', 'Conditional Death')