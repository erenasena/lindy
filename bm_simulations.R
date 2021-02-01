### Outline: 
# 1-Brownian motion simulations
# 2-Describing the hitting time distributions
# 3-Computing hazard functions and conditional survival (the Lindy effect) 

### The necessary libraries
library(dplyr)
library(parallel) # to run the parallel simulations 
library(ggplot2) # to plot the Brownian motion simulations
library(reshape) # for visualizations
library(survival) # survival analysis
library(survminer) # to plot the survival curves 
library(dynpred) # to compute conditional survival 
library(qualityTools) # for the Q-Q plots 
library(evir) # for the mean excess plot 
library(ineq) # for the Zenga plot 
library(fitdistrplus) # to fit distributions
library(logspline) # to test K-S hypothesis
library(actuar)
library(statmod) # to work with inverse Gaussian values 

### The Brownian motion functions 

## Arithmetic Brownian Motion (resamples until the value is below the reflecting barrier)  
my_abm <- function(nsim, t0, t, n, X0, mu, sigma, L, R){ # nsim is the number of sample paths 
                                                         # t0 is the first time point (e.g., 0)
                                                         # t is the final time point (e.g., 1)
                                                         # n is the number of time intervals that divide t0-t
                                                         # mu = drift, sigma = diffusion, L = absorb. barrier, R = reflecting barrier 
  
  dt <- t/n # the size of the time intervals, also the variance of BM increments 
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
        
        if(X[i,j] <= R){ # if we are equal or below, the value stays the same
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
    values <- list("values" = X, "time" = event_time, "event" = event_status)
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
    values <- list("values" = X, "time" = event_time, "event" = event_status)
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
  u <- u[-which(names(u) == "time")]
  u <- u[-which(names(u) == "event")]  
  m <- matrix(data = u, nrow = nsim, ncol = length(0:n), byrow = TRUE)
  df <- data.frame(x = t(m), row.names = paste0("time", 0:n, ""))
  colnames(x = df) <- paste0("sim", 1:nsim, "")
  store <- list(m, df)
  return(store)
}

# A function to store the hitting times 
times <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "time")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("sim", 1:nsim, ""))
  colnames(x = df) <- "time"
  store <- list(m, df)
  return(store)
}

# A function to store the event / censoring times 
events <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "event")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("sim", 1:nsim, ""))
  colnames(x = df) <- "event"
  store <- list(m, df)
  return(store)
}

## Parallel runs 
f <- function(i){ # specify the desired function and parameter values here
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 0, mu = -1, sigma = 0.5, L = -0.5, R = 10000000)
  #my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = -10000000, R = 100000000) 
}

set.seed(1)
res <- mclapply(X = 1:100, f, mc.cores = 8, mc.set.seed = TRUE)

v <- values(x = res, nsim = 100, n = 1000) # indexing the BM values 
m_val <- v[[1]] # BM values in a matrix (goes into the plotting function)
df_val <- v[[2]] # BM values in a data frame

t <- times(x = res, nsim = 100, n = 1000) # indexing the hitting times 
m_times <- t[[1]] # in a matrix (for histograms)
df_times <- t[[2]] # in a data frame 

e <- events(x = res, nsim = 100, n = 1000)
m_event <- e[[1]] 
df_event <- e[[2]]

data <- cbind(df_times, df_event)

### Plotting all types of BM 

## With no absorbing barrier  
bm <- m_val
abm <- m_val
gbm <- m_val

## Plots

# Some common params
time = c(1:1001)
cl <- rainbow(100)

# Standard BM 
pdf(file = "sbm.pdf")
plot(x = time, y = bm[1, ], type = "l", col = cl[1], 
     main = "Standard Brownian motion", xlab = "Time", ylab = "Values", 
     ylim = c(-3, 3), lty = 1, las = 1, bty = "n")
invisible(lapply(2:100, function(i) lines(x = time, y = bm[i, ], col = cl[i], type = 'l')))
dev.off()

# ABM
pdf(file = "abm.pdf")
plot(x = time, y = abm[1, ], type = "l", col = cl[1], 
     main = "Arithmetic Brownian motion", xlab = "Time", ylab = "Values", 
     ylim = c(-3, 1), lty = 1, las = 1, bty = "n")
invisible(lapply(2:100, function(i) lines(x = time, y = abm[i, ], col = cl[i], type = 'l')))
dev.off()

# GBM 
pdf(file = "gbm.pdf")
plot(x = time, y = gbm[1, ], type = "l", col = cl[1], 
     main = "Geometric Brownian motion", xlab = "Time", ylab = "Values", 
     ylim = c(0, max(gbm)), lty = 1, las = 1, bty = "n")
invisible(lapply(2:100, function(i) lines(x = time, y = gbm[i, ], col = cl[i ], type = 'l')))
dev.off()

## BM with an absorbing barrier 
bmb <- m_val
B <- -0.1 # barrier level 
pdf(file = "bmb.pdf")
plot(x = time, y = bmb[1, ], type = "l", col = cl[1], 
     main = "BM with an absorbing barrier", xlab = "Time", ylab = "Values", 
     ylim = c(-0.5, 2), lty = 1, las = 1, bty = "n")
invisible(lapply(2:100, function(i) lines(x = time, y = bmb[i, ], col = cl[i], type = 'l')))
abline(h = B, col = "red") # add the absorbing barrier
dev.off()

## ABM with an absorbing barrier 
babm <- m_val
B <- -0.5
pdf(file = "babm.pdf")
plot(x = time, y = babm[1, ], type = "l", col = cl[1], 
     main = "ABM with an absorbing barrier", xlab = "Time", ylab = "Values", 
     ylim = c(-0.7, 0.5), lty = 1, las = 1, bty = "n")
invisible(lapply(2:100, function(i) lines(x = time, y = babm[i, ], col = cl[i], type = 'l')))
abline(h = B, col = "red") # add the absorbing barrier
dev.off()

# Data 
#std <- readRDS(file = 'std') # standard Brownian motion 
#abm1 <- readRDS(file = 'abm1')
#gbm1 <- readRDS(file = 'gbm1') # gbm with abs but no ref 
#gbm2 <- readRDS(file = 'gbm2') # gbm with both barriers
#highnets <- readRDS(file = 'networks') # networks with high connectivity 
#lownets <- readRDS(file = 'lownets') # networks with low connectivity 

#data <- lownets
#length(which(data$time == 1)) / length(data$time)
#which(data$time > 70)
### Distribution fitting 

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

bmplot(x = abm, nsim = 1000, n = 10000, L = 100000, R = 100000, 
       ylim = c(0, 1000), title = "Geometric Brownian motion with an absorbing")

### Describing the distributions

# Some distributions
set.seed(12345)
g <- rgamma(n = 100000, shape = 0.5, rate = 1)
w <- rweibull(n = 100000, shape = 0.9, scale = 50)
e <- rexp(n = 100000, 1)
p <- rpareto(n = 100000, scale = 1, shape = 3)

# Histogram of the hitting times
hist(data$time, breaks = 10, xlim = c(0, 100000), main = 'Geometric Brownian motion 
     with absorbing and reflecting barriers', xlab = "Time", col = "lightblue", border = "darkblue", prob = F)

# Q-Q Plot to check exponentiality (if linear, thin tails, if concave, there may be heavy tailedness)
#hittings <- sort(data$time) # sort the data
hittings <- sort(w)
p <- ppoints(hittings, length(hittings)) # get the probabilities of the data
s <- quantile(x = hittings, p = p) # sample quantiles
q <- qexp(p = p) # exponential quantiles 
qqplot(x = s, y = q, xlab = "Sample quantiles", ylab = "Theoretical quantiles", 
       main = "Exponential QQ Plot")

# Zipf / log-log plot to check for power law decay (linearity indicates power law)
hittings <- sort(data$time) # sort the data
fit <- ecdf(hittings) # empirical cumulative distribution function
s <- 1 - fit(hittings) # empirical survival function  
logs <- log(s) # the log of the survival function 
logx <- log(hittings) # the log of the sorted failure times 

plot(x = logx, y = logs, xlab = "log(failure times)", ylab = "log(survival function)", 
     main = "The Zipf Plot of the Weibull distribution", bty = 'n', xlim = c(min(logx), max(logx)))
legend(x = "bottomleft", legend = c('c = 0.8:1.3', 'delta = 0.1', 'threshold = 5', 
                                't = 10000', 'nsim = 1000')) 

# The Zipf plot function of Cirillo 
zipfplot <- function (data, type = "plot", title = TRUE){
  # type should be equal to ’points’ if you want to add the
  # Zipf Plot to an existing graph
  # With other strings or no string a new graph is created.
  # If title is set to be F, the title of the plot is not given.
  # This can be useful when embedding the Zipf plot into other
  # plots.
  data <- sort(as.numeric(data)) #sorting data 
  y <- 1 - ppoints(data) #computing 1-F(x)
  if (type == "points") {
    points(data, y, xlog = T, ylog = T, xlab = "x on log scale",
           ylab = "1-F(x) on log scale")
    }
    plot(data, y, log = "xy", xlab = "x on log scale", ylab = "1-F(x) on log scale", 
         main = "Zipf Plot of the Pareto Distribution", bty = "n")
}

zipfplot(data = g)

png(filename = "zipfweibull")
zipfplot(data = w)
dev.off()

png(filename = "zipfpareto")
zipfplot(data = p)
dev.off()

### Survival analysis 

## Create survival objects and fit hazards for all Brownian motions 
surv_data <- data.frame(Time = data$time, Event = data$event, row.names = paste0("Sim", 1:length(data$time), ""))
surv_object <- Surv(time = data$time, event = data$event) 
surv_fit <- survfit(surv_object ~ 1, data = data, ctype = 1)
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)

survival <- surv_fit$surv
surv_time <- surv_fit$time

plot(x = haz_fit$time, y = haz_fit$hazard, xlab = "Time", 
     ylab = "Hazard", ylim = c(min(haz_fit$hazard), max(haz_fit$hazard)),
     type = 'l', bty = 'n', main = 'Hazard function of GBM stopping times', col = 'blue')

## Conditional survival
fit <- dynpred::Fwindow(object = surv_fit, width = 100, variance = TRUE, conf.level = 0.95)
condeath <- fit$Fw
consurv <- 1 - condeath
plot(x = fit$time, y = consurv, type = 'l', col = 'green', xlab = 'Time', 
     ylab = 'Conditional Probability', main = 'Conditional survival and death
     probabilities of GBM stopping times',
     ylim = c(0, 1), xlim = c(0, 1000), bty = 'n')
lines(x = fit$time, y = condeath, col = 'red') # change the limit of the y-axis to c(0, 1) to see this 
legend(x = 'bottom', bty = 'n', lty = c(1, 1), 
       col = c("green", "red"), 
       legend = c("Conditional survival", "Conditional death"))

#cond <- data.frame(con_time, con_surv, con_death)
#colnames(cond) <- c('Time', 'Conditional Survival', 'Conditional Death')
