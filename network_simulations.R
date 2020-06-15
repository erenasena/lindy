### Libraries
library(parallel)
library(IsingFit) # make sure you install version 0.3.0 (NOT 0.3.1) from the source:
#install.packages("/Users/erenasena/Downloads/IsingFit_0.3.0.tar", repos = NULL, type = "source")
sessionInfo() # check if the correct version is attached 

### Load the data 
md_data <- data.matrix(read.table("md_data.txt", header = F))
colnames(md_data) <- c("dep.mood", "loss.int", "w.loss", "w.gain", "dec.app", "inc.app", "insomnia", 
                       "hypsomnia", "psycm.agit", "psycm.ret", "fatigue", "worth", "conc.prob", "death")


################## The Functions ################## 

# This function computes the total number of active symptoms at each time point for n networks. 
lindynet <- function(data, c, t) { # data must be a matrix 
                                   # c is the connectivity parameter 
                                   # t is the number of time points

  fit <- IsingFit::IsingFit(x = data, plot = FALSE) # x must be cross-sectional data
  
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) 
  W <- fit$weiadj 
  P <- X <- matrix(data = 0, nrow = t, ncol = j) # matrix with probabilities (P) and simulated data (X): initially 0s. 
  W <- c * W 
  
  for (t in 2:t) { 
    A <- X[t - 1,] %*% W # the activation function (formula (1) in paper)
    P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active ( = 1) if the probability is greater 
                                    # than randomly chosen uniformly distributed number between 0 & 1
  }
  states <- apply(X, 1, sum) # compute the total number of symptoms per row (time point)
  network <- list("states" = states, "symptoms" = X)
  return(network)
}

# This function prepares the data for survival analysis 
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
  
################## The Simulations ################## 

## Set up parallel runs 
RNGkind("L'Ecuyer-CMRG") # this is the random number generator needed in parallel processing 
detectCores() # tells you the number of cores your computer can use for the simulations 

## Set up the parameters 
t <- 10000
nsim <- 100
c <- 1.1
threshold <- 5 # how many active symptoms to consider the network "depressed"? 
n <- 1 # how many successive time points does a symptom need to be 0 to be considered "dead"? 
j <- 14 # number of symptoms

## Prepare the function 
f <- function(i) {
  lindynet(data = md_data, c = c, t = t)
}

## Run the simulations and store the data 
set.seed(1)
results <- mclapply(X = 1:nsim, FUN = f, mc.cores = 8, mc.set.seed = TRUE)

# The network states 
states <- paste0("states", 1:t, "") # total number of active symptoms in each time point
states <- matrix(data = unlist(results)[which(names(unlist(results)) == states)], nrow = nsim, ncol = t, byrow = TRUE)
dimnames(x = states) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:t, ""))
data <- as.data.frame(survnet(X = states, threshold = threshold)) # transition times
hist(data$time, breaks = 15, xlab = 'Time', main = 'Time distribution of transitions')
legend(x = "center", legend = c(paste("c", "=", c), paste("time", "=", time), 
                                paste("nsim", "=", nsim), paste("threshold", "=", threshold)))

# The symptom states
symptom <- matrix(data = NA, nrow = t, ncol = ncol(md_data))
symptoms <- rep(list(x = symptom), nsim)
for(i in 1:nsim) {
  symptoms[[i]] <- t(results[[i]]$symptoms)
  dimnames(x = symptoms[[i]]) <- list( c("dep.mood", "loss.int", "w.loss", "w.gain", "dec.app", 
                                        "inc.app", "insomnia", "hypsomnia", "psycm.agit", 
                                        "psycm.ret", "fatigue", "worth", "conc.prob", "death"), 
                                       paste0("t", 1:t, ""))
}
names(x = symptoms) <- paste0("Network", 1:nsim, "")

# A function for the symptom level 
symptom <- function(L) {
  events <- matrix(nrow = nsim, ncol = 14)
  times <- matrix(nrow = nsim, ncol = 14)
  
  for(l in 1:length(L)) {
    X <- L[[l]]
    
    start <- numeric(length = nrow(X))
    for(i in 1:nrow(X)) {
      start[i] <- which(X[i,] == 1)[1] # index when the network is depressed for the first time
      }
    time <- numeric(length = nrow(X))
    event <- numeric(length = nrow(X))
    
    for(i in 1:nrow(X)) {
      for(j in start[i]:ncol(X)) { # start counting time from the first depressed point on 
        if(X[i, j] == 1 & j < ncol(X)) {
          time[i] <- time[i] + 1
          } else if(X[i, j] == 1 & j == ncol(X)) {
            time[i] <- time[i]
            event[i] <- 0
            } else if(X[i, j] == 0) {
              time[i] <- time[i]
              event[i] <- 1
              break 
            }
        }
    }
    events[l,] <- event
    dimnames(events) <- list(paste0("Network", 1:nsim, ""), 
                             c("dep.mood", "loss.int", "w.loss", "w.gain", "dec.app", "inc.app", 
                               "insomnia", "hypsomnia", "psycm.agit", "psycm.ret", "fatigue", 
                               "worth", "conc.prob", "death"))
    times[l,] <- time
    dimnames(times) <- list(paste0("Network", 1:nsim, ""), 
                                 c("dep.mood", "loss.int", "w.loss", "w.gain", "dec.app", "inc.app", 
                                   "insomnia", "hypsomnia", "psycm.agit", "psycm.ret", "fatigue", 
                                   "worth", "conc.prob", "death"))
  }
  data <- list('events' = events, 'times' = times)
  return(data)
}

data <- symptom(L = symptoms)

dep_times <- cbind(data$times[,'dep.mood'])
dep_event <- cbind(data$events[,'dep.mood'])

lossint_times <- cbind(data$times[, 'loss.int'])
lossint_event <- cbind(data$events[, 'loss.int'])

decapp_times <- cbind(data$times[,'dec.app'])
decapp_events <- cbind(data$events[,'dec.app'])

w.loss_times <- cbind(data$times[,'w.loss'])
w.loss_events <- cbind(data$events[,'w.loss'])

fatigue_times <- cbind(data$times[,'fatigue'])
fatigue_events <- cbind(data$events[,'fatigue'])

worth_times <- cbind(data$times[,'worth'])
worth_events <- cbind(data$events[,'worth'])

death_times <- cbind(data$times[,'death'])
death_events <- cbind(data$events[,'death'])

conc_times <- cbind(data$times[,'conc.prob'])
conc_events <- cbind(data$events[,'conc.prob'])

insom_times <- cbind(data$times[,'insomnia'])
insom_events <- cbind(data$events[,'insomnia'])

hist(dep_times, breaks = 25, main = "Failure time distribution of loss of interest", xlab = "Time")
legend(x = "center", legend = c(paste("c", "=", c), paste("time", "=", t), paste("nsim", "=", nsim)))
data <- data.frame('time' = dep_times, 'event' = dep_event)

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
