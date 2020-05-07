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

### The BM functions 

## Arithmetic Brownian Motion (with drift and absorbing barrier) 
my_abm <- function(nsim, t0, t, n, X0, mu, sigma, L){ # same as above, except now we have mu for drift 
  # sigma is the diffusion coefficient / volatility 
  # also added a t0 term, which is just 0 in practice
  
  dt <- t/n 
  time <- seq(from = t0, to = t, by = dt) # a time vector that determines the number of columns, with rows being the simulated sample paths  
  initial <- X0  
  X <- matrix(nrow = nsim, ncol = length(time)) 
  X[1:nsim, 1] <- X0
  event_time <- numeric(length = nsim) # I changed the hitting times to event time and 
  event_status <- numeric(length = nsim) # added an event status indicating if failure happened or not
  
  for(i in 1:nrow(X)){
    for(j in 2:length(time)){
      X[i,j] <- X0 + mu * dt + sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1) 
      if(X[i,j] > L & j < ncol(X)){ 
        X0 <- X[i,j] 
      } else if(X[i,j] > L & j == ncol(X)){ 
        X0 <- initial 
        event_time[i] <- ncol(X) 
        event_status[i] <- 0 # censored, did not die until the last point 
      } else if(X[i,j] <= L){ 
        X0 <- initial 
        event_time[i] <- j 
        event_status[i] <- 1 # hit the boundary, dead 
        break
      }
    }
    values <- list("Values" = X, "Event time" = event_time, "Event status" = event_status)
  }
  return(values)
}

## Geometric Brownian Motion (with absorbing barrier)
# The only difference between this and the ABM function is the formula. The variables are the same. 
my_gbm <- function(nsim, t0, t, n, X0, mu, sigma, L){ 
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
bmplot <- function(x, nsim, n, L, ylim, title){ # x is the matrix output of the BM functions, n is the number of simulations, t is the vector of time points, as in the BM functions 
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
    ggtitle(title) + xlab("Time") + ylab("Value") + ylim(ylim) + geom_hline(yintercept = L, color = "orange", size = 0.3)
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
  my_gbm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90) 
  #my_abm(nsim = 1, t0 = 0, t = 1, n = 1000, X0 = 100, mu = -1, sigma = 1, L = 90)
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

#p <- bmplot(x = m_val, nsim = 1000, n = 1000, L = 90, ylim = c(87.5, 100), # Define the range of the y-axis  
            #title = "Brownian motion with an absorbing barrier")
#print(p)

# Histogram of hitting times
hist(m_times, breaks = c(seq(from = 0, to = 1010, by = 10)), xlim = c(0, 1000))

### Getting the hazard functions 
surv_data <- data.frame(Time = m_times, Event = m_event, row.names = paste0("Sim", 1:nrow(m_times), ""))
surv_object <- Surv(time = m_times, event = m_event) 
fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)
hazard <- plot(fit$time, fit$hazard, xlab='Time', ylab = 'Hazard Rate', type = 'l', xlim = c(0, 1000), ylim = c(min(fit$haz), max(fit$haz)))

### Hazard function trials

# The hazard function of Pareto
library(VGAM)
hazard <- function(x, a, b){
  fx <- dpareto(x, shape = a, scale = b)
  dx <- ppareto(x, shape = a, scale = b)
  sx <- 1-dx
  hx <- fx/sx
  return(hx)
}
plot(hazard(x = 1:100, a = 2, b = 1), type = "l")

# The hazard function of exponential 
eh <- function(x, b){
  for(i in 1:100)
    x[i] <- 1/b
  return(x)
}
evalues <- eh(x = 1:100, b = 1)
plot(evalues, type = "l")

# The hazard function of Weibull 
wh <- function(a, b, t){
  h <- (b/a)*((t/a)^(b-1))
  return(h)
}
wvalues <- wh(a = 1, b = 2, t = 1:100)
plot(wvalues, type = 'l')

## bshazard  
library(bshazard)

# Pareto: should be decreasing and it is
set.seed(1234)
prt <- VGAM::rpareto(n = 1000, shape = 2, scale = 1)
fit_1 <- bshazard(Surv(prt) ~ 1)
plot(x = fit_1$time, y = fit_1$hazard, ylim =c(min(fit_1$haz), max(fit_1$haz)), 
     xlab = "Time", ylab = "Hazard", type = 'l')

# Exponential: should be constant but it is not
set.seed(1234)
exps <- rexp(n = 100, 1)
fit_2 <- bshazard(Surv(exps) ~ 1)
plot(x = fit_2$time, y = fit_2$hazard, ylim = c(min(fit_2$haz), max(fit_2$haz)), 
     xlab = "Time", ylab = "Hazard", type = 'l')

# Weibull: should be decreasing for shape < 1, constant for shape = 1, increasing for shape > 1
# except for constant hazard, it's all good 
set.seed(1234)
wei <- rweibull(n = 1000, shape = 4, scale = 1)
fit_3 <- bshazard(Surv(wei) ~ 1)
plot(x = fit_3$time, fit_3$hazard, ylim = c(min(fit_3$haz), max(fit_3$haz)), 
     xlab = "Time", ylab = "Hazard", type = 'l')
