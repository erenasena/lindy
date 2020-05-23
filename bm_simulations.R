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

### Dynamic prediction

## Survival curve and fitting the model 
surv_data <- data.frame(Time = m_times, Event = m_event, row.names = paste0("Sim", 1:nrow(m_times), ""))
surv_object <- Surv(time = m_times, event = m_event) 
surv_fit <- survfit(surv_object ~ 1)
surv_plot <- plot(x = surv_fit$time, y = surv_fit$surv, type = 'l', xlab = "Time")

## Conditional death with the package 
con_death <- Fwindow(object = surv_fit, width = 1, variance = T, conf.level = 0.95)
deaths <- con_death$Fw # conditional deaths 
con_plot <- plot(x = con_death$time, y = deaths, type = 'l', ylim = c(min(con_death$Fw), max(con_death$Fw)),
     xlab = "Time", ylab = "Conditional death")

fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data) # fit the hazard function 
hazard <- plot(fit$time, fit$hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 1010), ylim = c(min(fit$haz), max(fit$haz)))

## My conditional survival function with the formula 
con_surv <- function(x){
  cons <- numeric(length = length(x) - 1)
  for(i in 1:length(cons)){
    cons[i] <- x[i+1] / x[i]
  }
  return(cons)
}

survs <- surv_fit$surv # survival probabilities 
surv_prob <- con_surv(x = survs) # conditional survival probabilities 
plot(x = 1:length(con), y = con, type = 'l', xlab = 'Time', ylab = 'Conditional survival')

## A function that computes the difference between successive rates 
diff <- function(x){
  diff <- numeric(length = length(x) - 1)
  for(i in 2:length(x)){
    diff[i-1] <- x[i] - x[i-1]
  }
  return(diff)
}

stat <- stationary_probs(x = hazard)
plot(x = 1:length(stat), y = stat, type = 'l', xlab = 'Index', ylab = 'Difference')
mean(stat)
median(stat)

### Generate perfect data from the exponential, Pareto and Weibull distributions 
set.seed(1)
pareto <- VGAM::rpareto(n = 10000, scale = 1, shape = 1.5)

set.seed(1)
exps <- rexp(n = 10000, rate = 1)

set.seed(1)
weibull <- rweibull(n = 1000, shape = 2)


surv_object <- Surv(time = exps) 
surv_fit <- survfit(surv_object ~ 1)
survs <- surv_fit$surv
con <- con_surv(x = survs)
plot(x = 1:length(con), y = con, type = 'l', xlab = 'Time', ylab = 'Conditional survival')

fit <- bshazard::bshazard(surv_object ~ 1)
hazard <- fit$hazard
hazard_plot <- plot(fit$time, fit$hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 50), ylim = c(min(fit$haz), max(fit$haz)))

