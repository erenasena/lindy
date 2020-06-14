### Libraries
library(IsingFit)
library(parallel)

### Load the data 
md_data <- data.matrix(read.table("md_data.txt", header = F))
colnames(md_data) <- c("dep", "int", "los", "gai", "dap", "iap", "iso", "hso", "agi", 
                       "ret", "fat", "wor", "con", "dea")

### This function computes the total number of active symptoms at each time point for n networks. 
lindynet <- function(data, c, n, nsim) { # data must be a matrix 
                                         # c is the connectivity parameter 
                                         # n is the number of time points
  
  states <- function(data, c, n) {

  fit <- IsingFit::IsingFit(x = data) # x must be cross-sectional data
  
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) 
  W <- fit$weiadj 
  P <- X <- matrix(0, n, j) # matrix with probabilities (P) and simulated data (X): initially 0s. 
  W <- c * W 
  
  for (t in 2:n) { 
    A <- X[t - 1,] %*% W # the activation function (formula (1) in paper)
    P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active ( = 1) if the probability is greater 
                                    # than randomly chosen uniformly distributed number between 0 & 1
  }
  states <- apply(X, 1, sum) # compute the total number of symptoms per row (time point)
  return(states)
  }
  
  f <- function(i) { # replicating the state function for n networks
    states(data = data, c = c, n = n)
  }
  res <- matrix(data = replicate(n = nsim, expr = states(data = data, c = c, n = n)), 
                  nrow = nsim, ncol = n, byrow = TRUE) 
  dimnames(res) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:n, ""))
  return(res)
}

set.seed(1)
X <- lindynet(data = md_data, c = 1.3, n = 10, nsim = 2) 

### This function prepares the data for survival analysis 
survnet <- function(X, threshold) { # X is the output matrix of the lindynet function 
                                    # threshold is the required number of active symptoms
  start <- numeric(length = nrow(X))
  for(i in 1:nrow(X)) {
    start[i] <- which(X[i,] >= threshold)[1] # index when the network is depressed for the first time
  }
  
  if(any(is.na(start)) == TRUE) { 
    X <- X[-(which(is.na(start) == TRUE)), ] # remove the networks that never became depressed
  } else { 
    X <- X 
  }
  
  start <- start[!is.na(start)]
  time <- numeric(length = nrow(X))
  event <- numeric(length = nrow(X))
  
  for(i in 1:nrow(X)) {
    for(j in start[i]:ncol(X)) { # start counting time from the first depressed point on 
      
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

data <- as.data.frame(survnet(X = X, threshold = 1))
hist(data$time, breaks = 10, xlab = 'Time', main = "Time distribution of transitions")
legend(x = "center", legend = c('c = 1.3', 'n = ', 'nsim = ', 'threshold = '))

## A function for the symptom level 
symptom <- function(X, n){
  event <- numeric(length = ncol(X))
  dead_start <- numeric(length = ncol(X))
  dead_end <- numeric(length = ncol(X))
  data_from_birth <- numeric()
  rle_data <- numeric()
  values <- numeric() # Values in the symptom: 0-1-0
  reps <- numeric() # How many times each value repeats 
  start <- numeric(length = ncol(X))
  
  for(j in 1:ncol(X)){ # First indexing the jth column of X 
    data <- cbind(X[,j])
    start[j] <- start_fun(data) # When the symptom is 1 for the first time in column j of X 
  }
  for(j in 1:length(start)){
    data_from_birth <- data[start[j]:nrow(data),]
    rle_data <- rle(as.vector(data_from_birth))
    values <- rle_data[[2]] 
    reps <- rle_data[[1]]
  }
  for(i in 1:length(values)){
    for(j in 1:length(event)){
      if(values[i] == 0 & reps[i] >= 2){
        event[j] <- 1
        break
      }
    }
  }
  for(i in 1:length(values)){
    for(j in 1:length(event)){
      if(values[i] == 0 & !(any(reps[i] >= 2))){
        event[j] <- 0
      }
    }
  }
  
  rep_deaths <- which(values == 0 & reps >= 2)
  deaths.lengths.cumsum <- cumsum(reps)
  ends <- deaths.lengths.cumsum[rep_deaths] # Until which time point the symptom remained 0. 
  newindex <- ifelse(rep_deaths > 1, rep_deaths-1, 0)
  dead <- deaths.lengths.cumsum[newindex] + 1
  if (0 %in% newindex){
    dead <- c(1, dead)} #The start point of death (0).
  
  for(j in 1:length(start)){
    if(event == 1){
      time <- start[j] + dead - 1 # The point at which it starts being 0 in succession 
      time_2 <- start[j] + ends - 1 # The last 0 in the succession of 0s
    } else {
      time <- nrow(s1)
      time_2 <- nrow(s1)
    }
  }
  return(dead_end)
}

symptom <- function()








































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
