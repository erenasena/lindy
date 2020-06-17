### Libraries
library(parallel)
library(IsingFit) # make sure you install version 0.3.0 (NOT 0.3.1) from the source:
#install.packages("/Users/erenasena/Downloads/IsingFit_0.3.0.tar", repos = NULL, type = "source")
#sessionInfo() # check if the correct version is attached 

### Load the data 
md_data <- data.matrix(read.table("md_data.txt", header = F))
colnames(md_data) <- c("dep.mood", "loss.int", "w.loss", "w.gain", "dec.app", "inc.app", "insomnia", 
                       "hypsomnia", "psycm.agit", "psycm.ret", "fatigue", "worth", "conc.prob", "death")

################## The Functions ################## 
# This function computes the total number of active symptoms at each time point for n networks. 
lindynet <- function(data, t, s1, s2, ds, c1, c2, dc) { # data must be a matrix 
                                                        # t is the number of time points
                                                        # s1 is the minimum stress
                                                        # s2 is the maximum stress
                                                        # ds is the difference between successive stress values
                                                        # c1 is the minimum connectivity
                                                        # c2 is the maximum connectivity 
                                                        # dc is the difference between successive connectivity values
  
  # Vectors to be sampled from in each network 
  stress <- seq(from = s1, to = s2, by = ds)
  connectivity <- seq(from = c1, to = c2, by = dc) # smaller increments are generally better
  
  # The Ising model  
  fit <- IsingFit::IsingFit(x = data, plot = FALSE) # x must be cross-sectional data
  
  # The network
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) 
  W <- fit$weiadj 
  P <- X <- matrix(data = 0, nrow = t, ncol = j) # matrix with probabilities (P) and simulated data (X): initially 0s. 
  c <- sample(x = connectivity, size = 1, replace = FALSE) # the connectivity parameter is random for each network
  W <- c * W 

  for (t in 2:t) { 
    if(t == 2) { # if the network is at the randomly sampled time point, 
      s <- sample(x = stress, size = 1, replace = FALSE) # a randomly sampled stress parameter
    } else {
      s <- 0 
    }
    A <- X[t - 1,] %*% W + s 
    P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active ( = 1) if the probability is greater 
                                    # than randomly chosen uniformly distributed number between 0 & 1
  }
  states <- apply(X, 1, sum) # compute the total number of symptoms per row (time point)
  network <- list("states" = states, "symptoms" = X)
  return(network)
}

# This function prepares the data for survival analysis 
survnet <- function(X, threshold) { # X is 'states', the output matrix of the lindynet function 
                                    # threshold is the required number of active symptoms
  time <- numeric(length = nrow(X))
  event <- numeric(length = nrow(X))
  
  for(i in 1:nrow(X)) { # for each network
    for(j in 2:ncol(X)) { # start counting time when the networks are depressed 
      if(X[i, j] >= threshold & j == ncol(X)) {
        event[i] <- 0
        time[i] <- j - 1

      } else if(X[i, j] < threshold) {
        event[i] <- 1
        time[i] <- j - 1
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
data <- md_data
nsim <- 1000
threshold <- 5 # this is arbitrary but a reasonable point 
t <- 10000 # 0-2000 time points ALWAYS show Lindy; actually 20%, even if everything has failed by 1000 points 
s1 <- 1 # 1 - 5; the lower end gives more barely depressed networks, 5 starts everything at 14 
s2 <- 15 # not sure about this one 
ds <- 1 # 1 is good 
c1 <- 0.90 # 1.15 and 1.1 gives good results; smaller and we don't have enough medium values, everything fails
c2 <- 1.05 # 1.25 - 1.30; smaller and we don't have enough large values, bigger and nothing fails 
dc <- 0.00015 # smaller values are better 

## Prepare the function 
f <- function(i) {
  lindynet(data = data, t = t, s1 = s1, s2 = s2, ds = ds, c1 = c1, c2 = c2, dc = dc)
}

## Run the simulations and store the data 
set.seed(1)
results <- mclapply(X = 1:nsim, FUN = f, mc.cores = 8, mc.set.seed = TRUE)

# The network states 
states <- paste0("states", 1:t, "") # total number of active symptoms in each time point
states <- matrix(data = unlist(results)[which(names(unlist(results)) == states)], nrow = nsim, ncol = t, byrow = TRUE)
dimnames(x = states) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:t, ""))
#states <- states[which(states[,2] >= threshold),] # include only the networks that became depressed when hit by stress
nrow(states) / nsim # the proportion of networks kept 

# Create the data frame for survival analysis
data <- as.data.frame(survnet(X = states, threshold = threshold)) # transition times
hist(data$time, breaks = 25, xlab = 'Time', main = 'Time distribution of transitions')
legend(x = "right", legend = c('1000 networks', '10,000 time points', 
                                  'threshold = 5', 'Min. stress = 1', 'Max. stress = 15',
                                  'delta stress = 1', 'Min. conn. = 0.90', 'Max. conn. = 1.05', 
                                  'delta conn. = 0.00015'), bty = 'n')

### Miscellaneous functions

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
