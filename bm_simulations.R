### Outline: 
# 1-Brownian motion simulations
# 2-Describing the hitting time distributions
# 3-Computing hazard functions and conditional survival (the Lindy effect) 

### The necessary libraries
library(dplyr) # don't know why 
library(parallel) # to run the parallel simulations 
library(ggplot2) # to plot the Brownian motion simulations
library(survival) # survival analysis
library(survminer) # to plot the survival curves 
library(dynpred) # to compute conditional survival 
library(qualityTools) # for the Q-Q plots 
library(evir) # for the mean excess plot 
library(ineq) # for the Zenga plot 

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
  my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90, R = 11000000000) 
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 0, mu = 0, sigma = 1, L = -0.1, R = 10000000)
}

set.seed(1)
res <- mclapply(X = 1:100, f, mc.cores = 8, mc.set.seed = TRUE)

v <- values(x = res, nsim = 100, n = 1000) # indexing the BM values 
m_val <- v[[1]] # BM values in a matrix (goes into the plotting function)
df_val <- v[[2]] # BM values in a data frame

t <- times(x = res, nsim = 100, n = 1000) # indexing the hitting times 
m_times <- t[[1]] # in a matrix (for histograms)
df_times <- t[[2]] # in a data frame 

e <- events(x = res, nsim = 100, n = 10000)
m_event <- e[[1]] 
df_event <- e[[2]]

data <- cbind(df_times, df_event)

### Describing the hitting time distribution

# Maximum to sum (MS) plot; converges to 0 if the moments are defined. Does not make sense to get 
# moments if they are infinite. So checking this first. 
ms <- function(x, p){ # x is the hitting times vector, p is the moment for which you want to test 
  ms <- numeric(length = length(x))
  n <- numeric(length = length(x))
  for(i in 1:length(x)){
    max <- max(x[1:i] ^ p)
    sum <- sum(x[1:i] ^ p)
    ms[i] <- max / sum
    n[i] <- i
  }
  return(list("MS" = ms, "n" = n))
}

ratio <- ms(x = data$time, p = 4) # if kurtosis is defined, the rest is defined also 
plot(x = ratio$n, y = ratio$MS, type = 'l', xlab = "Number of values", ylab = "Max / sum ratio")

# Descriptive statistics of the hitting times; moments are defined so we can get mean and sd 
mean <- mean(data$time)
sd <- sd(data$time)
quantile(data$time) 
ext <- data$time[which(data$time > 500)] # the extreme values 
length(ext) / length(data$time) # what proportion of data are larger than a certain value 
sum(ext) / sum(data$time) # what proportion of the sum they make

# Histogram of the hitting times
hist(data$time, breaks = 100, xlim = c(0, 10000), main = 'GBM with an absorbing barrier')
legend(x = "center", legend = c('c = 0.8:1.3', 'delta = 0.1', 'threshold = 5', 
                                't = 10000', 'nsim = 1000'))

# Q-Q Plot to check exponentiality (if linear, thin tails, if concave, there may be heavy tailedness)
hittings <- sort(data$time) # sort the data
p <- ppoints(hittings, length(hittings)) # get the probabilities of the data
s <- quantile(x = hittings, p = p) # sample quantiles
q <- qexp(p = p) # exponential quantiles 
qqplot(x = s, y = q, xlab = "Sample quantiles", ylab = "Theoretical quantiles", 
       main = "Exponential QQ Plot")

qqPlot(x = data$time, y = "exponential", xlab = "Sample quantiles", 
       ylab = "Theoretical quantiles", main = "Exponential Q-Q Plot") # shorter way with confidence bands

# Zipf / log-log plot to check for power law decay (linearity indicates power law)
hittings <- sort(data$time) # sort the data
fit <- ecdf(hittings) # empirical cumulative distribution function
s <- 1 - fit(hittings) # empirical survival function  
logs <- log(s) # the log of the survival function 
logx <- log(hittings) # the log of the sorted failure times 

plot(x = logx, y = logs, xlab = "log(failure times)", ylab = "log(survival function)", main = "The Zipf Plot") 
legend(x = "bottomleft", legend = c('c = 0.8:1.3', 'delta = 0.1', 'threshold = 5', 
                                't = 10000', 'nsim = 1000'))                              
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
    points(data, y, xlog=T, ylog=T, xlab = "x on log scale",
           ylab = "1-F(x) on log scale")
    }
  else if(title==F){
    plot(data, y, log="xy", xlab = "x on log scale", ylab = "1-F(x) on log scale")
    }
  else{
    plot(data, y, log="xy", xlab = "x on log scale", ylab = "1-F(x) on log scale", main="Zipf Plot")
    }
}

zipfplot(data = data$time)

# Mean excess (ME) plot (linearity indicates power law, concavity lognorm, constant exp, decreasing norm)
evir::meplot(sort(data$time)) 
VGAM::meplot(sort(data$time)) # gives confidence bands
legend(x = "bottomleft", legend = c('c = 1.225', 'n = 10000', 'nsim = 1000'))

meplot <- function(data, cut = 5){
  # In cut you can specify the number of maxima you want to exclude. # The standard value is 5
  data <- sort(as.numeric(data)); n = length(data);
  mex <- c();
  for (i in 1:n){
    mex[i] <- mean(data[data > data[i]]) - data[i];
  }
  data_out <- data[1:(n - cut)];
  mex_out <- mex[1:(n - cut)];
  plot(data_out, mex_out, xlab = "Threshold u", ylab = "Mean Excess e(u)", main = "Mean Excess Plot (Meplot)")
}

meplot(data = data$time, cut = 5)

# Discriminant moment ratio plot; this is supposed to discriminate the distribution but is a bit faulty
moment_plot <- function(data){
  # "data" is a vector containing the sample data
  ############################################## ############################################## 
  # CV and Skewness functions 
  coefvar <- function(data){
    CV <- sd(data)/mean(data)
    CV
  }
  skewness <- function(data){
    m_3 <- mean((data - mean(data)) ^ 3) 
    skew <- m_3 / (sd(data) ^ 3)
    skew
  }
############################################## ############################################## 
  # Computation of CV and Skewness
  # CV
  CV <- coefvar(data); 
  # Skewness 
  skew <- skewness(data) 
  # Rule of Thumb
  if (CV < 0 | skew < 0.15){print("Possibly neither nor lognormal. Thin tails."); stop}
############################################## # Preparation of the plot ############################################## 
  ############################################## # Paretian Area
  # The upper limit - Pareto I 
  p <- seq(3.001, 400, length.out = 250) 
  g2brup <- 1 / (sqrt(p * (p - 2))) 
  g3brup <- (1 + p) / (p - 3) * 2 / (sqrt(1 - 2 / p))
  # The lower limit, corresponding to the Inverted Gamma 
  g2ibup <- seq(0.001, 0.999, length.out = 250) 
  g3ibup <- 4 * g2ibup / (1 - g2ibup ^ 2) 
  ##############################################
  # Lognormal area
  # Upper limit: Lognormal
  w <- seq(1.01, 20, length.out = 250)
  g2log <- sqrt(w - 1)
  g3log <- (w + 2) * sqrt(w - 1)
  # Lower limit - Gamma
  g2iblow <- seq(0, 20, length.out = 250)
  g3iblow <- 2 * g2iblow
  ##############################################
  # Exponential Area
  # The upper limit corresponds to the lower limit of the
  # lognormal area
  # The lower limit - Bernoulli
  g2below <- seq(0, 20, length.out = 250)
  g3below <- g2below - 1 / g2below 
  # The Gray area is obtained for free from
  # the previous lines of code. 
  # Normal / Symmetric distribution
  g2nor <- seq(0, 20, length.out = 250)
  g3nor <- rep(0, 250)

  # PLOT
  # Limits 
  plot(g2iblow, g3iblow, "l", xlab = "CV", ylab = "Skewness", main = "Discriminant Moment-ratio Plot", xlim = c(0, 20), ylim = c(-1, 40)) 
  lines(g2ibup, g3ibup, "l")
  lines(g2brup, g3brup, "l")
  lines(g2below, g3below, "l")
  lines(g2log, g3log, lty = 2) # Lognormal
  lines(g2nor, g3nor, lty = 2) # Normal
  # Strictly Paretian Area 
  polygon(c(g2ibup, g2brup), c(g3ibup, g3brup), col = "green") 
  points(0, 2, pch = 1, cex = 0.8) # Pareto limit point
  # Hints for interpretation
  text(-0.2, 20, cex = 0.8, srt = 90, "Pareto I") 
  text(1.2, 20, cex = 0.8, srt = 90, "Inverted Gamma")
  text(2.5, 12, cex = 0.8, srt = 70, "Lognormal") 
  text(12, 21, cex = 0.8, srt = 23, "Gamma") 
  text(14, 11, cex = 0.8, srt = 10, "Bernoulli") 
  text(15, 1.5, cex = 0.8, "Normal or Symmetric") 
  points(CV, skew, pch = 16, col = "red")
  points(CV, skew, pch = 16, col = "red")
  return(c(CV, skew))
}

moments <- moment_plot(data = data$time)
legend(x = "top", legend = c('c = 1.225', 'n = 10000', 'nsim = 1000'))

# Zenga plot 
zengaplot <- function(data){
  # Since the code relies on the Lorenz curve
  # as computed by the "ineq" library,
  # we upload it
  library(ineq)
  # Empirical Lorenz
  est <- Lc(data)
  # Zenga curve
  Zu <- (est$p - est$L) / (est$p * (1 - est$L))
  # We rescale the first and the last point for
  # graphical reasons
  Zu[1] <- Zu[2]; Zu[length(Zu)] <- Zu[(length(Zu)-1)]
  # Here’s the plot
  plot(est$p, Zu, xlab = "u", ylab = "Z(u)", ylim = c(0, 1), main = "Zenga plot", "l", lty = 1)
}

zengaplot(data = data$time)
legend(x = "center", legend = c('c = c(0.85, 1.3, 1.25)', 'delta = 0.01', 'threshold = 5', 
                                't = 10000', 'nsim = 1000'))

# Fit a distribution with various methods
fitdistrplus::fitdist(data = data$time, distr = "lnorm", method = "mle") # parameter est. 
ks.test(x = data$time, y = "plnorm", 'meanlog' = 1.733794, 'sdlog' = 2.346450)
lnorm <- rlnorm(n = 100, 1.733794, 2.346450) # trying to see if matches the other graphs

### Survival analysis 
surv_data <- data.frame(Time = data$time, Event = data$event, row.names = paste0("Sim", 1:length(data$time), ""))
surv_object <- Surv(time = data$time, event = data$event) 
surv_fit <- survfit(surv_object ~ 1)
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)

survival <- surv_fit$surv
surv_time <- surv_fit$time

hazard <- haz_fit$hazard
haz_time <- haz_fit$time

hazard_plot <- plot(x = haz_time, y = hazard, xlab = 'Time', ylab = 'Hazard Rate', type = 'l', 
                    xlim = c(min(haz_time), max(haz_time)), ylim = c(min(haz_fit$haz), max(haz_fit$haz)))

## Conditional survival
fit <- dynpred::Fwindow(object = surv_fit, width = 100, variance = TRUE, conf.level = 0.95)
con_time <- fit$time # the calculated times
con_death <- fit$Fw # conditional death 
con_surv <- 1 - con_death # conditional survival; they are mirror images
plot(x = con_time, y = con_surv, type = 'l', col = 'green', xlab = 'Time', 
     ylab = 'Probability', main = 'Conditional survival and death over time',
     ylim = c(0, 1), xlim = c(min(con_time), max(con_time)))
lines(con_time, con_death, col = 'red') # change the limit of the y-axis to c(0, 1) to see this 
cond <- data.frame(con_time, con_surv, con_death)
colnames(cond) <- c('Time', 'Conditional Survival', 'Conditional Death')
