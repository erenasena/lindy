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
lindynet <- function(data, t, C, c1, c2, dc, s1, s2, ds) { # data must be a matrix 
                                                           # t is the number of time points
                                                           # If C = T, c is sampled from a vector for each network; else it's c1 for all networks
                                                           # c1 is the minimum connectivity
                                                           # c2 is the maximum connectivity 
                                                           # dc is the difference between successive connectivity values
                                                           # s1 is the minimum stress
                                                           # s2 is the maximum stress
                                                           # ds is the difference between successive stress values
  
  # The Ising model  
  fit <- IsingFit::IsingFit(x = data, plot = FALSE) # x must be cross-sectional data
  
  # The network
  j <- ncol(data) # number of symptoms
  b <- abs(fit$thresholds) 
  W <- fit$weiadj 
  P <- X <- matrix(data = 0, nrow = t, ncol = j) # matrix with probabilities (P) and simulated data (X): initially 0s. 
  
  if(C == TRUE) {
    c <- sample(x = seq(from = c1, to = c2, by = dc), size = 1, replace = FALSE)
    } else {
      c <- c1
      }
  
  W <- c * W 
  
  for (t in 2:t) {
    if(t == 2) { # if the network is at the randomly sampled time point, 
      s <- sample(x = seq(from = s1, to = s2, by = ds), size = 1, replace = FALSE) # a randomly sampled stress parameter
      } else { 
        s <- 0 
      }
    
    A <- X[t - 1,] %*% W + s 
    P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
    X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active ( = 1) if the probability is greater 
                                    # than randomly chosen uniformly distributed number between 0 & 1
  }
  states <- apply(X, 1, sum) # compute the total number of symptoms per row (time point)
  network <- list("states" = states, "connectivity" = c, "symptoms" = X)
  return(network)
}

# This function prepares the data for survival analysis 
survnet <- function(X, threshold) { # X is 'states', the output matrix of the lindynet function 
                                    # threshold is the required number of active symptoms
                                    # connectivity is the connectivity vector, output from lindynet
  time <- numeric(length = nrow(X))
  event <- numeric(length = nrow(X))
  
  for(i in 1:nrow(X)) { # for each network
    for(j in 2:ncol(X)) { # start counting time after hitting with stress 
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
t <- 10000 # 0-2000 time points ALWAYS show Lindy; actually 20%, even if everything has failed by 1000 points 

C <- T
c1 <- 1.00 # 1.15 and 1.1 gives good results; smaller and we don't have enough medium values, everything fails
c2 <- 1.15 # 1.25 - 1.30; smaller and we don't have enough large values, bigger and nothing fails 
dc <- 0.00015 # smaller values are better 

s1 <- 1 # 1 - 5; the lower end gives more barely depressed networks, 5 starts everything at 14 
s2 <- 15 # not sure about this one 
ds <- 1 # 1 is good 

threshold <- 5 # this is arbitrary but a reasonable point 

## Prepare the function 
f <- function(i) {
  lindynet(data = data, t = t, C = C, c1 = c1, c2 = c2, dc = dc, s1 = s1, s2 = s2, ds = ds)
}

## Run the simulations and store the data 
set.seed(1)
results <- mclapply(X = 1:nsim, FUN = f, mc.cores = 8, mc.set.seed = TRUE)

# The network states 
states <- paste0("states", 1:t, "") # total number of active symptoms in each time point
states <- matrix(data = unlist(results)[which(names(unlist(results)) %in% states)], 
                 nrow = nsim, ncol = t, byrow = TRUE)
dimnames(x = states) <- list(paste0("Network", 1:nsim, ""), paste0("t", 1:t, ""))
#states <- states[which(states[,2] >= threshold),] # include only the networks that became depressed when hit by stress
nrow(states) / nsim # the proportion of networks kept 

# The network connectivity levels 
connectivity <- unlist(results)[which(names(unlist(results)) == "connectivity")]
