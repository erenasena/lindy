### Libraries 
library(IsingFit)
library(qgraph)

### Load the data 
data <- read.table("md_data.txt", header = F)
md_data <- data.matrix(data)
colnames(md_data) <- c("dep", "int", "los", "gai", "dap", "iap", "iso", "hso", "agi", 
                       "ret", "fat", "wor", "con", "dea")
### Simulations

## Creating a function that returns the states over time for each connectivity level
state <- function(data, c, n){ # data is a matrix 
                               # c is the connectivity vector to be looped over
                               # n is the number of time points
  
  ################ Creating the network ################ 
  
  library(IsingFit)
  library(qgraph)
  fit <- IsingFit(md_data) # x must be cross-sectional data
  pdf(file = "IsingModelKendlerData.pdf")
  qgraph(fit$weiadj, layout = "spring", filetype = "pdf", filename = "IsingModelKendlerData")
  
  ################ Set up the visuals ################ 
  
  pdf(file = "ResultsBasicModel.pdf", width = 10, height = 2.5)
  par(mfrow = c(1,3), col.main = "white")
  
  ################ Creating the network ################ 
  
    j <- ncol(data) # number of symptoms
    b <- abs(fit$thresholds) # thresholds 
    W <- fit$weiadj # weights 
    P <- X <- matrix(0, n, j) # matrix with probabilities (P) and simulated data (X); full of 0s initially
    W <- c * W # now W has 3 matrices, each with the original weights multiplied by a connectivity value
    
    for (t in 2:n) { 
      A <- X[t-1,] %*% W # the activation function (formula (1) in paper)
      P[t,] <- 1 / (1 + exp(b - A)) # the probability function (formula (2) in paper)
      X[t,] <- 1 * (P[t,] > runif(j)) # symptom i becomes active (=1) if the probability is greater than randomly chosen uniformly distributed number between 0 and 1
    }
  state <- apply(X, 1, sum)
  return(state)
}
  
replicate(2, state(data = md_data, c = 1.3, n = 100))

## Some plotting stuff
#if (c == 1) {StateWeak <- State}
#plot(State[1:1500], type="l", xlim = c(1, 1500), ylim = c(0, 14), axes = F, ann = F, col = "grey")
#axis(side = 1, c(1, 500, 1000, 1500))
#if (c == 1) {axis(side = 2, c(0, 7, 14))} else {axis(side = 2, at = c(0, 7, 14), labels = F)}


