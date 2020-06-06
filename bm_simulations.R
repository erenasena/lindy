### Survival and hazard functions of hitting times 

# The necessary libraries
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
  my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90, R = 1100000) 
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = 0, sigma = 0.15, L = 9000000, R = 10000000)
}

set.seed(1)
res <- mclapply(X = 1:1000, f, mc.cores = 8, mc.set.seed = TRUE)

v <- values(x = res, nsim = 1000, n = 1000) # indexing the BM values 
m_val <- v[[1]] # BM values in a matrix (goes into the plotting function)
df_val <- v[[2]] # BM values in a data frame

t <- times(x = res, nsim = 1000, n = 1000) # indexing the hitting times 
m_times <- t[[1]] # in a matrix (for histograms)
df_times <- t[[2]] # in a data frame 

e <- events(x = res, nsim = 1000, n = 1000)
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
time <- haz_fit$time
hazard_plot <- plot(time, hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', 
                    xlim = c(0, 1000), ylim = c(min(haz_fit$haz), max(haz_fit$haz)))



## Polynomial regression

# Prepare the data 
data <- data.frame(time = time, hazards = hazard)
timelims <- range(time)
time.grid <- seq(from = timelims[1], to = timelims[2])

# The models 
fit.1 <- lm(hazard ~ time, data = data)
fit.2 <- lm(hazard ~ poly(time, 2), data = data)
fit.3 <- lm(hazard ~ poly(time, 3), data = data)
fit.4 <- lm(hazard ~ poly(time, 4), data = data)
fit.5 <- lm(hazard ~ poly(time, 5), data = data)
fit.6 <- lm(hazard ~ poly(time, 6), data = data)
fit.7 <- lm(hazard ~ poly(time, 7), data = data)
fit.8 <- lm(hazard ~ poly(time, 8), data = data)
fit.9 <- lm(hazard ~ poly(time, 9), data = data)
fit.10 <- lm(hazard ~ poly(time, 10), data = data)
fit.20 <- lm(hazard ~ poly(time, 20), data = data)

# Plotting the predictions
preds <- predict(fit.10, newdata = list(time = time.grid), se = TRUE)
se.bands <- cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)

plot(x = time, y = hazard, xlim = timelims, cex = .5, col = "darkgrey")
lines(time.grid , preds$fit, lwd = 2, col = "blue")
matlines(time.grid, se.bands , lwd = 1, col = "blue", lty = 3)

# Comparing the models 
anova(fit.1, fit.2, fit.3, fit.4, fit.5, fit.6, fit.7, fit.8, fit.9, fit.10, fit.20)

## Regression splines
library(splines)

# Natural spline 
fit <- lm(hazard ~ ns(time, df = 8), data = data)
plot(x = time, y = hazard, xlim = timelims, cex = .5, col = "darkgrey")
pred = predict(fit, newdata = list(time = time.grid), se = T)
lines(time.grid, pred$fit, col= "black", lwd = 2)
lines(time.grid, pred$fit + 2 * pred$se.fit, lty = "dashed")
lines(time.grid, pred$fit - 2 * pred$se.fit, lty = "dashed")
summary(fit)

# Smoothing splines 
fit2 <- smooth.spline(time, hazard, cv = TRUE)
plot(time, hazard, xlim = timelims, cex = .5, col = "darkgrey")
lines(fit2, col = "red", lwd = 2)

# Comparing models 


