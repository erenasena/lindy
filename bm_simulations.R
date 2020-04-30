#### Brownian Motion Simulations

# The section "Simulations", which is at the bottom of the file, has the parallel runs and two functions to store values. 
# Make sure you run these two functions as well as the functions "my_bm", "my_gbm", and "bmplot" before the parallel simulations. 

### Install the required packages
#install.packages("fExpressCertificates")
#install.packages("somebm")
#install.packages("yuima")

### Load the function packages
library(fExpressCertificates)
library(somebm)
library(yuima)

### My functions 

## Standard Brownian Motion

# This function is actually not needed for our project because it does not have a drift. Also, 
# the ABM function does the same thing with mu = 0, provided we set the absorbing barrier L 
# sufficiently above or below 0. This is because in the standard BM, the starting value is always 0. 
# Consequently, the process may hit L - despite not having drift - and set the remaining values to NA. 
# This will not be a problem in ABM and GBM, because we will start with a value above 0). 

my_bm <- function(t, n, nsim){ # t is the terminal/final time point (e.g., the end of 1 year)
                               # n is the number of intervals between 0:t (e.g., 365 days)
                               # nsim is the number of sample paths (i.e., different BM sims)
  sig2 <- t/n # the variance of each increment equals its length
  time <- seq(from = 0, to = t, by = sig2) # time vector with t/sig2 + 1 = n + 1 elements 
  sigma <- sqrt(sig2) # the standard deviation of each increment 
  ri <- rnorm(n = nsim * (length(time) - 1), mean = 0, sd = sigma) # simulating 100 normally distributed increments for all simulations
  X <- matrix(data = ri, nrow = nsim, ncol = length(time) - 1) # storing these increments in a i(number of simulations) x j(number of increments in a simulation) matrix 
  X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum))) # combining nsim 0s for the starting value of each BM and
  return(X)                                        # the new values of the BM as the cumulative sum of the value of the BM at each time point
}

## Arithmetic Brownian Motion (with drift and absorbing barrier) 
my_abm <- function(nsim, t, n, X0, mu, sigma, L){ # nsim is the number of BM sample paths to be simulated 
                                                  # t is the final time point 
                                                  # n is the number of time intervals/increments from 0:t
                                                  # X0 is the first value of BM at time 0, mu is the drift 
                                                  # sigma is the diffusion coefficient and L is the barrier
  
  dt <- t/n # time divided into equal intervals; change in time and also the variance of the increments
  sig2 <- sigma^2 # the volatility term by which the process is scaled, not the variance of the increments 
  time <- seq(from = 0, to = t, by = dt) # a time vector that determines the number of columns, with rows being the simulated sample paths  
  
  initial <- X0 # the initial value of each simulated sample path 
  X <- matrix(nrow = nsim, ncol = length(time)) # the matrix that will store all the values
  X[1:nsim, 1] <- X0 # setting the first value of each simulation 
  hit_times <- numeric(length = nsim)
  
  for(i in 1:nrow(X)){
    for(j in 2:length(time)){
      X[i,j] <- X0 + mu * dt + sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1) # the formula for the next value of BM
      if(X[i,j] > L & j < ncol(X)){ # if the value is larger than the absorbing barrier L, meaning we have not hit it yet, and if the time point (ncol) is not the last time point of a given simulation,
        X0 <- X[i,j] # then the initial value of BM is updated by the new value 
      } else if(X[i,j] > L & j == ncol(X)){ # if we have not hit L but we have come to the last time point in the simulation,
        X0 <- initial # then the initial value is updated by 1, which is the first value of every simulation
        hit_times[i] <- length(time) # and the hit time is none, so I set it as the last time point
      } else if(X[i,j] <= L){ # if we have hit the absorbing barrier from above in a given simulation
        X0 <- initial # then the initial value is updated by 1, since we will start the new BM, and
        hit_times[i] <- j # the time at which we reach L (the column number) updates the hitting time vector
        break
      }
    }
    values <- list(Values = X, Hittings = hit_times)
  }
  return(values)
}

## Geometric Brownian Motion (with absorbing barrier)
# The only difference between this and the ABM function is the formula. The variables are the same. 
my_gbm <- function(nsim, t, n, X0, mu, sigma, L){ # nsim is the number of BM sample paths to be simulated 
                                                  # t is the final time point 
                                                  # n is the number of time intervals/increments from 0:t
                                                  # X0 is the first value of BM at time 0, mu is the drift 
                                                  # sigma is the diffusion coefficient, L is the barrier
  
  dt <- t/n # time divided into equal intervals; change in time and also the variance of the increments
  sig2 <- sigma^2 # the volatility term by which the process is scaled, not the variance of the increments 
  time <- seq(from = 0, to = t, by = dt) # a time vector that determines the number of columns, with rows being the simulated sample paths  
  
  initial <- X0 # the initial value of each simulated sample path 
  X <- matrix(nrow = nsim, ncol = length(time)) # the matrix that will store all the values
  X[1:nsim, 1] <- X0 # setting the first value of each simulation 
  hit_times <- numeric(length = nsim) # one hitting time for each simulated sample path  
  
  for(i in 1:nrow(X)){
    for(j in 2:length(time)){
      X[i,j] <- X0 * exp(((mu - 0.5 * sig2) * dt) + (sigma * sqrt(dt) * rnorm(n = 1, mean = 0, sd = 1))) # the formula for the next value of BM
      if(X[i,j] > L & j < ncol(X)){ # if the value is larger than the absorbing barrier L, meaning we have not hit it yet, and if the time point (ncol) is not the last time point of a given simulation,
        X0 <- X[i,j] # then the initial value of BM is updated by the new value 
      } else if(X[i,j] > L & j == ncol(X)){ # if we have not hit L but we have come to the last time point in the simulation,
        X0 <- initial # then the initial value is updated by 1, which is the first value of every simulation
        hit_times[i] <- length(time) # and the hit time is none
      } else if(X[i,j] <= L){ # if we have hit the absorbing barrier from above in a given simulation
        X0 <- initial # then the initial value is updated by 1, since we will start the new BM, and
        hit_times[i] <- j # the time at which we reach L (the column number) updates the hitting time vector
        break
      }
    }
    values <- list(Values = X, Hittings = hit_times)
  }
  return(values)
}

## GBM with the yuima package. This mostly gives the same results, except when "Terminal"/the last point
# is much larger (e.g., 100 instead of 1 for all functions). This is not a problem with the other packages. 
require(yuima)
set.seed(1)
m <- setModel(drift = "mu*s", diffusion = "sigma*s", state.var = "s", time.var = "t", solve.var = "s", xinit = 1)
options(warn = -1)# Suppress Warnings from R
simnum <- 100 # number of BMs to be simulated 
grid <- setSampling(Initial = 0, Terminal = 1, n = 100) # Initial is the first time point, Terminal is the last time point, n is the number of intervals 
newsim <- function(i){simulate(m, sampling = grid, true.param = list(mu = -1, sigma = 1))@data@original.data}
options(warn = -1)# Suppress Warnings from R
sim <- sapply(1:simnum, function(x)newsim(x))
gbm4 <- t(sim) # this is the final matrix with the values of all the simulations 

# This is independent from the package. I added it to get the hitting times
#yuima_hittings <- numeric(length = nrow(gbm4)) 
#for(i in 1:nrow(gbm4)){
  #for(j in 2:ncol(gbm4)){
    #if(gbm4[i,j] <= 0.001){
      #yuima_hittings[i] <- j
      #break
    #}
  #}
#}

### Comparisons

## Standard Brownian motion
set.seed(1)
bm1 <- my_bm(t = 1, n = 100, nsim = 100) 
set.seed(1)
bm2 <- my_abm(nsim = 100, t = 1, n = 100, X0 = 0, mu = 0, sigma = 1, L = -100) # my abm with mu = 0
set.seed(1)
bm3 <- fExpressCertificates::BM(S0 = 0, mu = 0, sigma = 1, T = 1, N = 100) # other abm with mu = 0
set.seed(1)
bm4 <- somebm::bm(x0 = 0, t0 = 0, t = 1, n = 100)

# Indexing the first 10 values for comparison 
SBM <- data.frame(MyBM = bm1[1, 1:10], MyABM = bm2[[1]][1, 1:10], fExpress = bm3[1:10], somebm = bm4[1:10])

## Arithmetic Brownian motion with drift
set.seed(1)
abm1 <- my_abm(nsim = 100, t = 1, n = 100, X0 = 1, mu = -1, sigma = 1, L = -100)
set.seed(1)
abm2 <- fExpressCertificates::BM(S0 = 1, mu = -1, sigma = 1, T = 1, N = 100)

# Indexing the first 10 values for comparison 
ABM <- data.frame(MyABM = abm1[[1]][1, 1:10], fExpress = abm2[1:10])

## Geometric Brownian motion: All except the yuima package give the same result (and yuima is very similar)
set.seed(1)
gbm1 <- my_gbm(nsim = 100, t = 1, n = 100, X0 = 1, mu = -1, sigma = 1, L = -100) # my gbm 
set.seed(1)
gbm2 <- fExpressCertificates::GBM(S0 = 1, mu = -1, sigma = 1, T = 1, N = 100)
set.seed(1)
gbm3 <- somebm::gbm(x0 = 1, mu = -1, sigma = 1, t0 = 0, t = 1, n = 100)
gbm4 <- gbm4 # the yuime output; this is defined above because the yuima package has its own thing (not important)

# Indexing the first 10 values for comparison 
GBM <- data.frame(MyGBM = gbm1[[1]][1, 1:10], fExpress = gbm2[1:10], somebm = gbm3[1:10], yuima = gbm4[1, 1:10])

## Comparing if the values are identical when rounded up to 8 decimals 
# apply(round(GBM, 8), 2, identical, round(GBM[,1], 8))
# apply(round(ABM, 8), 2, identical, round(ABM[,1], 8))
# apply(round(SBM, 8), 2, identical, round(SBM[,1], 8))

## Printing the data frames defined above 
#print(SBM)
#print(ABM)
#print(GBM)

### Visualizations
library(ggplot2)
library(reshape2)
bmplot <- function(x, nsim, t, n, L, ylim, title){ # x is the matrix output of the BM functions, n is the number of simulations, t is the vector of time points, as in the BM functions 
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


#bmplot(x = abm1[[1]], nsim = 100, t = 1, n = 100, ylim = c(min(abm1[[1]]), max(abm1[[1]])), title = "Arithmetic Brownian motion with negative drift and an absorbing barrier")
#bmplot(x = bm1, nsim = 100, t = 1, n = 100, ylim = c(min(bm1), max(bm1)), title = "Arithmetic Brownian motion with negative drift and an absorbing barrier")
#bmplot(x = bm2[[1]], nsim = 100, t = 1, n = 100, ylim = c(min(bm2[[1]]), max(bm2[[1]])), title = "Arithmetic Brownian motion with negative drift and an absorbing barrier")

### Simulations
library(parallel)
RNGkind("L'Ecuyer-CMRG") # this is the random number generator needed in parallel processing 
detectCores() # tells you the number of cores your computer can use for the simulations 

## Storing the values: These two functions will allow us to nealty store the resutls. It is easier to 
## define them before running the simulations. I wrote them because the results of parallel simulations are a mess. 

# The BM values 
values <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) == "Hittings")]
  m <- matrix(data = u, nrow = nsim, ncol = length(0:n), byrow = TRUE)
  df <- data.frame(x = t(m), row.names = paste0("Time", 0:n, ""))
  colnames(x = df) <- paste0("Sim", 1:nsim, "")
  store <- list(m, df)
  return(store)
}

# The hitting times 
hittings <- function(x, nsim, n){
  u <- unlist(x)
  u <- u[-which(names(u) != "Hittings")]
  m <- matrix(data = u, nrow = nsim, ncol = 1, byrow = TRUE)
  df <- data.frame(x = m, row.names = paste0("Sim", 1:nsim, ""))
  colnames(x = df) <- "Hitting time"
  store <- list(m, df)
  return(store)
}

## Parallel runs 

run_parallel <- function(args, function_choice){
  if(function_choice == "abm"){ # specify the desired funciton here 
    f <- function(i){ # the parameter values here
      my_abm(nsim = args$nsim, t = args$t, n = args$n, X0 = args$xo, mu = args$mu, sigma = args$sigma, L = args$L)
    }
} else {
    f <- function(i){
      my_gbm(nsim = args$nsim, t = args$t, n = args$n, X0 = args$xo, mu = args$mu, sigma = args$sigma, L = args$L)
    }
  }
}

set.seed(1)
res <- mclapply(X = 1:100, f, mc.cores = 8, mc.set.seed = TRUE) # X is the n of sim as a vector 
                                                                # f is the function defined above 
                                                                # mc.cores is the number of cores you want to use 

s <- values(x = res, nsim = 100, n = 100) # indexing the BM values 
mval <- s[[1]] # BM values in a matrix (goes into the plotting function)
dfval <- s[[2]] # BM values in a data frame

h <- hittings(x = res, nsim = 100, n = 100) # indexing the hitting times 
mhit <- h[[1]] # in a matrix (for histograms)
dfhit <- h[[2]] # in a data frame 

p <- bmplot(x = mval, nsim = 100, t = 1, n = 100, L = 0, ylim = c(-0.5, 1), # Define the range of the y-axis  
             title = "Arithmetic Brownian motion with negative drift and an absorbing barrier")
print(p)

hist(mhit) # Histogram of hitting times 
    
}
