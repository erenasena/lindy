### Libraries 
library(IsingFit)
library(parallel)

### Load the data 
md_data <- data.matrix(read.table("md_data.txt", header = F))
colnames(md_data) <- c("dep", "int", "los", "gai", "dap", "iap", "iso", "hso", "agi", 
                       "ret", "fat", "wor", "con", "dea")
### Simulations

## Creating a function that returns the states over time for each connectivity level.
## No stress parameter. 
state <- function(data, c, n) { # data is a matrix 
                                      # c is the connectivity parameter 
                                      # n is the number of time points
  
  ################ Fitting the Ising model to estimate weights and biases ################ 
  
  fit <- IsingFit(data = data) # x must be cross-sectional data
  
  ################ Creating the network ################ 
  
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
  state <- apply(X, 1, sum)
  return(state)
}

state(data = md_data, c = 1.2, n = 10)

# Running the simulations in parallel 
lindynet <- function(data, c, n, nsim) { # connectivity, number of time points, number of networks to be simulated
  
  RNGkind("L'Ecuyer-CMRG") # this is the random number generator needed in parallel processing 

  f <- function(i) { 
    state(data = data, c = c, n = n)
  }
  
  res <- unlist(mclapply(1:nsim, FUN = f, mc.cores = 8, mc.set.seed = TRUE))
  states <- matrix(data = res, nrow = nsim, ncol = n, byrow = TRUE) 
  dimnames(states) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:n, ""))
  return(states)
}

set.seed(1)
states <- lindynet(data = md_data, c = 1.2, n = 4, nsim = 2)



## This function adds the stress parameter
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