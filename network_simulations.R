### Libraries
library(IsingFit)
library(parallel)

### Load the data 
md_data <- data.matrix(read.table("md_data.txt", header = F))
colnames(md_data) <- c("dep", "int", "los", "gai", "dap", "iap", "iso", "hso", "agi", 
                       "ret", "fat", "wor", "con", "dea")
### Simulations
states <- function(data, c, n) { # data is a matrix 
  # c is the connectivity parameter 
  # n is the number of time points
  
  ####### Fitting the Ising model to estimate weights and biases ####### 
  
  fit <- IsingFit::IsingFit(x = data) # x must be cross-sectional data

    ####### Creating the network ####### 
  
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) # thresholds 
  W <- fit$weiadj # weights 
  P <- X <- matrix(0, n, j) # matrix with probabilities (P) and simulated data (X); full of 0s initially
  W <- c * W # multiply the weights by the connectivity parameter 
  
  for (t in 2:n) { 
    A <- X[t - 1,] %*% W # the activation function (formula (1) in paper)
    P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active ( = 1) if the probability is greater than randomly chosen uniformly distributed number between 0 and 1
  }
  states <- apply(X, 1, sum)
  return(states)
}
## Creating a function that returns the number of symptoms per time point for n networks 
lindynet <- function(data, c, n, nsim) { # connectivity, number of time points, number of networks to be simulated
  
  
  ################ This function creates a network ################ 
  f <- function(i){
    states(data = data, c = c, n = n)
  }
  
  ################ This function replicates the state function and stores the values ################ 
  
  res <- matrix(data = replicate(n = nsim, expr = states(data = data, c = c, n = n)), 
                nrow = nsim, ncol = n, byrow = TRUE) 
  dimnames(res) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:n, ""))
  return(res)
}

set.seed(1)
X <- lindynet(data = md_data, c = 1.2, n = 10000, nsim = 100) 

# For each network, index when they are depressed for the first time 
start <- function(X, threshold) { # X is the output matrix of the lindynet function 
                                  # threshold is the required number of active symptoms for depression 
  start <- numeric(length = nrow(X))
  for(i in 1:nrow(X)) {
    start[i] <- which(X[i,] >= threshold)[1]
  }
  return(start)
}

# For each network, start counting time when it's on and stop when it goes below the threshold
index <- function(X, threshold) {
  start <- start(X = X, threshold = threshold)
  if(any(is.na(start)) == TRUE) {
    X <- X[-(which(is.na(start) == TRUE)), ] # remove the networks that never became depressed
    } else { 
      X <- X 
  }
  start <- start[!is.na(start)]
  time <- numeric(length = nrow(X))
  event <- numeric(length = nrow(X))
  
  for(i in 1:nrow(X)) {
    for(j in start[i]:ncol(X)) {
      
      if(X[i, j] >= threshold & j < ncol(X)) {
        time[i] <- time[i] + 1
        
      } else if(X[i, j] >= threshold & j == ncol(X)) {
        time[i] <- time[i]
        event[i] <- 0
        
      } else if(X[i, j] < threshold) {
        time[i] <- time[i]
        event[i] <- 1
        break 
      }
    }
  }
  data <- cbind(event, time)
  row.names(data) <- paste0("Network", 1:nrow(X), "")
  return(data)
}

failures <- index(X = X, threshold = 5)
data <- as.data.frame(failures)
hist(data$time)

### Survival analysis 
surv_data <- data.frame(Time = data$time, Event = data$event, row.names = paste0("Sim", 1:nrow(failures), ""))
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
fit <- dynpred::Fwindow(object = surv_fit, width = 1000, variance = TRUE, conf.level = 0.95)
con_time <- fit$time # the calculated times
con_death <- fit$Fw # conditional death 
con_surv <- 1 - con_death # conditional survival; they are mirror images
plot(x = con_time, y = con_surv, type = 'l', col = 'green', xlab = 'Time', 
     ylab = 'Probability', main = 'Conditional survival and death over time',
     ylim = c(0, 1), xlim = c(min(con_time), max(con_time)))
lines(con_time, con_death, col = 'red') # change the limit of the y-axis to c(0, 1) to see this 
cond <- data.frame(con_time, con_surv, con_death)
colnames(cond) <- c('Time', 'Conditional Survival', 'Conditional Death')

## This function adds the stress parameter but is incomplete
stress <- function(data, c, n){
  
  ################ Fitting the Ising model to estimate weights and biases ################ 
  fit <- IsingFit(data) # x must be cross-sectional data
    
  ################ Creating the network ################ 
  S <- delta <- rep(0, n)
  delta[1] <- 1
  step <- 0.01
  Smax <- 15
    
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) # thresholds 
  W <- fit$weiadj # weights 
  P <- X <- matrix(0, n, j) # matrix with probabilities (P) and simulated data (X); full of 0s initially
  W <- c * W # now W has 3 matrices, each with the original weights multiplied by a connectivity value
    
  for (t in 2:n) {
    if (abs(S[t - 1]) > Smax) {
      delta[t] <- -1 * delta[t - 1]
    } else {
        delta[t] <- delta[t - 1]
    }
    
    S[t] <- S[t - 1] + delta[t] * step
    A <- X[t - 1,] %*% W + S[t] #the activation function (formula (1) in paper)
    P[t,] <- 1 / (1 + exp(b - A)) #the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) #symptom i becomes active (=1) if the probability is greater than randomly chosen uniformly distributed number between 0 and 1
  }
  X <- X[-1:-100,]
  S <- S[-1:-100]
  delta <- delta[-1:-100]
  
  S <- 1 * S
  s <- apply(X, 1, sum)
  Sseq <- seq(-Smax, Smax, .2)
  Sc <- cut(S, Sseq)
  test <- tapply(s, list(Sc, delta), mean)
  #plot(Sseq[26:100], test[c(26:100), 1], type = "b", col = 'black', axes = F, ann = F)
  #lines(Sseq[26:100], test[26:100, 2], type = "b", col = 'grey', ann = F)
  
  if(c == 2) {
    mtext("Stress", side = 1, line = 3, col = "black", adj = 0.5)  
    axis(side = 1, c(-10, -5, 0, 5))
  }
  if (c == 1) {
    axis(side = 2, c(0, 7, 14))
  } else {
      axis(side = 2, at = c(0, 7, 14), labels = F)
  }
  
  if(c == 1) {
    sdata <- cbind(S, s) 
    } else {
      sdata <- cbind(sdata, s)
    }
  
  if(c == 1) {
    mdata <- cbind(Sseq[-1], tapply(s, list(Sc, delta), mean))
    } else {
      mdata <- cbind(mdata, tapply(s, list(Sc, delta), mean))
  }
  return()
}

stress(data = md_data, c = 1.2, n = 10)
