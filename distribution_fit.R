### Libraries 
library(philentropy) # Information theoretical methods 
library(fitdistrplus) # distribution fit 
library(statmod) # to work with the Inverse-Gaussian
library(VGAM) # to work with Levy 
library(latex2exp)

### Generating example hazard functions to display 
t <- c(1:100) # a time vector for all 

# The hazard function of Weibull 
w <- wh(a = 50, b = 4.5, t = t)

# The hazard function of exponential 
eh <- function(x, b){
  for(i in 1:100)
    x[i] <- 1/b
  return(x)
}
e <- eh(x = t, b = 2)

# Pareto hazard function from Eliazar (2017)
ph <- function(p, t){
  h <- ((1+p)/p)*(1/t)
  return(h)
}
p <- ph(p = 1/4, t = t)

# Plotting all the curves together 
pdf(file = "diff_hazards.pdf")
plot(x = t, y = w, type = "l", col = "red", lty = 1, ylim = c(0, 1), 
     xlab = "Time", ylab = "Hazard rate", main = "Hazard functions", las = 1, bty = "n")
lines(x = t, y = e, col = "blue", lty = 1, ylim = c(0, 1))
lines(x = t, y = p, col = "dark green", lty = 1, ylim = c(0, 1))
dev.off()

########### The Levy Distribution (hitting time of standard BM) ###########
n = 10000
#a = 0.1 # the distance between X0 and B 
#c = a^2

## Plotting the Levy density 

# Values to plot over 
x <- seq(from = 0.01, to = 3, length = n)

# Compute the densities 
d <- VGAM::dlevy(x = x, location = 0, scale = 1/2)

# Plot for different scale values
pdf(file = "levydensities.pdf")

plot(x, d, type = "l", lty = 1, xlab = "Values",
     ylab = "Density", main = "The Levy Distribution", 
     bty = "n", ylim = c(0, 1), col = "red", las = 1)

legend("topright", legend = c("scale = 1/2", "scale = 1", "scale = 2"),
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

c <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(c)){
  d <- VGAM::dlevy(x = x, location = 0, scale = c[i])
  lines(x, d, lwd = 1, col = col[i])
}

dev.off()

### Hazard function of Levy
levyhaz <- function(x, location, scale){ 
  h <- numeric(length = length(x))
  x <- sort(x)
  for(i in 1:length(h)){
    h[i] <- VGAM::dlevy(x = x[i], location, scale) / (1 - VGAM::plevy(q = x[i], location, scale))
  }
  return(h)
}

## Compute the hazard rates 

# Generate random values 
set.seed(12345)
levy <- VGAM::rlevy(n, location = 0, scale = 1/2)

# Compute hazards 
h <- levyhaz(x = levy, location = 0, scale = 1/2)

# Plot all for different scales 
pdf(file = "levyhazards.pdf")

plot(x = c(1:length(h)), y = h, type = 'l', 
     main = "The Hazard Function of the Levy Distribution", xlab = "Time", 
     ylab = "Hazard rate", ylim = c(0, max(h)), 
     bty = "n", col = "red", las = 1)

legend("topright", legend = c("scale = 1/2", "scale = 1", "scale = 2"),
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

c <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(c)){
  levy <- VGAM::rlevy(n, location = 0, scale = c[i])
  h <- levyhaz(x = levy, location = 0, scale = c[i])  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()
########### Working with the Inverse Gaussian (IG) ########### 

### Density of the Inverse Gaussian 

# Values to plot over 
x <- seq(from = 0, to = 3, length = 100000)

# Compute the densities 
d <- statmod::dinvgauss(x = x, mean = 1, shape = 1)

# Plot for different scale values
pdf(file = "IGdensities.pdf")

plot(x, d, type = "l", lty = 1, xlab = "Values",
     ylab = "Density", main = "The Inverse Gaussian Distribution", 
     bty = "n", ylim = c(0, max(d)), col = "red", las = 1)

legend("topright", legend = c(TeX("$\\mu = 1$, $\\lambda = 1$"), 
                              TeX("$\\mu = 1$, $\\lambda = 2$"), 
                              TeX("$\\mu = 1$, $\\lambda = 4$"), 
                              TeX("$\\mu = 2$, $\\lambda = 1$"), 
                              TeX("$\\mu = 4$, $\\lambda = 1$")), 
       col = c("red", "blue", "dark green", "black", "orange"), lty = 1, cex = 0.8, bty = "n")

mu <- c(1, 1, 2, 4)
lambda <- c(2, 4, 1, 1)
col <- c("blue", "dark green", "black", "orange")

for (i in 1:length(mu)){
  d <- statmod::dinvgauss(x = x, mean = mu[i], shape = lambda[i])
  lines(x, d, lwd = 1, col = col[i])
}

dev.off()

### The hazard function of the IG (general)

## The hazard function 
IGhaz <- function(x, mu, lambda){
  x <- sort(x, decreasing = FALSE)
  h <- numeric(length = length(x))
  for(i in 1:length(h)){
    h[i] <- statmod::dinvgauss(x = x[i], mean = mu, shape = lambda) / (1 - statmod::pinvgauss(q = x[i], mean = mu, shape = lambda))
  }
  return(h)
}

## Plotting the hazard function

# Generate random values 
set.seed(12345)
IG <- statmod::rinvgauss(n, m = 1, s = 1/2)
h <- IGhaz(x = IG, mu = 1, lambda = 1/2)

# Plot the function 
pdf(file = "IGhazards.pdf")

plot(x = c(1:length(h)), y = h, type = 'l', 
     main = "The Hazard Function of the Inverse Gaussian", 
     xlab = "Time", ylab = "Hazard rate", ylim = c(0, 2), 
     bty = "n", col = "red", las = 1)

legend("topright", legend = c(TeX("$\\mu = 1$, $\\lambda = 1/2$"), 
                              TeX("$\\mu = 1$, $\\lambda = 1$"), 
                              TeX("$\\mu = 2$, $\\lambda = 1$"), 
                              TeX("$\\mu = 4$, $\\lambda = 1$")), 
       col = c("red", "blue", "dark green", "orange"), lty = 1, cex = 0.8, bty = "n")

mu <- c(1, 2, 4)
lambda <- c(1, 1, 1)
col <- c("blue", "dark green", "orange")

for (i in 1:length(mu)){
  IG <- statmod::rinvgauss(n, mean = mu[i], shape = lambda[i])
  h <- IGhaz(x = IG, mu = mu[i], lambda = lambda[i])  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()

########### Lognormal ########### (decreasing for mu > 1, sigma > 1)

### Lognormal density 

# Values to plot over 
x <- seq(from = 0, to = 3, length = 100000)

# Compute the densities 
d <- dlnorm(x = x, meanlog = 0, sdlog = 1/2)

# Plot for different scale values
pdf(file = "lnormdensities.pdf")

plot(x, d, type = "l", lty = 1, xlab = "Values",
     ylab = "Density", main = "The Lognormal Distribution", 
     bty = "n", ylim = c(0, 1.5), col = "red", las = 1)

legend("topright", legend = c(TeX("$\\mu = 0$, $\\sigma = 1/2$"), 
                              TeX("$\\mu = 0$, $\\sigma = 1$"), 
                              TeX("$\\mu = 0$, $\\sigma = 2$")), 
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

sigma <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(sigma)){
  d <- dlnorm(x = x, meanlog = 0, sdlog = sigma[i])
  lines(x, d, lwd = 1, col = col[i])
}

dev.off()

### Lognormal hazard (abs value affected by mu but the shape only by sigma)
lnorm_haz <- function(x, mu, sigma){
  x <- sort(x, decreasing = FALSE)
  h <- numeric(length = length(x))
  for(i in 1:length(h)){
    h[i] <- dlnorm(x = x[i], meanlog = mu, sdlog = sigma) / (1 - plnorm(q = x[i], meanlog = mu, sdlog = sigma))
  }
  return(h)
}

# Generate random values 
set.seed(12345)
lnorm <- rlnorm(n = 10000, meanlog = 0, sdlog = 1/2)
h <- lnorm_haz(x = lnorm, mu = 0, sigma = 1/2)

# Plot the function 
pdf(file = "lnormhazards.pdf")

plot(x = c(1:length(h)), y = h, type = 'l', 
     main = "The Hazard Function of the Lognormal", 
     xlab = "Time", ylab = "Hazard rate", ylim = c(0, 2), 
     bty = "n", col = "red", las = 1)

legend("right", legend = c(TeX("$\\mu = 0$, $\\sigma = 1/2$"), 
                              TeX("$\\mu = 0$, $\\sigma = 1$"), 
                              TeX("$\\mu = 0$, $\\sigma = 2$")), 
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

sigma <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(sigma)){
  lnorm <- rlnorm(n = 10000, meanlog = 0, sdlog = sigma[i])
  h <- lnorm_haz(x = lnorm, mu = 0, sigma = sigma[i])  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()

########### Weibull ########### (decreasing for shape < 1)

### Weibull density 

# Values to plot over 
x <- seq(from = 0, to = 3, length = 100000)

# Compute the densities 
d <- dweibull(x = x, shape = 1/2, scale = 1)

# Plot for different scale values
pdf(file = "wdens.pdf")

plot(x, d, type = "l", lty = 1, xlab = "Values",
     ylab = "Density", main = "The Weibull Distribution", 
     bty = "n", ylim = c(0, 2), col = "red", las = 1)

legend("topright", legend = c(TeX("$\\lambda = 1$, $\\kappa = 1/2$"), 
                              TeX("$\\lambda = 1$, $\\kappa = 1$"), 
                              TeX("$\\lambda = 1$, $\\kappa = 3/2$"), 
                              TeX("$\\lambda = 1$, $\\kappa = 4$")), 
       col = c("red", "blue", "dark green", "orange"), lty = 1, cex = 0.8, bty = "n")

shape <- c(1, 3/2, 4)
col <- c("blue", "dark green", "orange")

for (i in 1:length(shape)){
  d <- dweibull(x = x, shape = shape[i], scale = 1)
  lines(x, d, lwd = 1, col = col[i])
}

dev.off()

### Weibull hazard 

# From the analytical form 
wh <- function(x, lambda, kappa){
  x <- sort(x, decreasing = FALSE)
  h <- (kappa/lambda)*((x/lambda)^(kappa-1))
  return(h)
}

# The same way as all the above 
weibull_haz <- function(x, shape, scale){
  x <- sort(x, decreasing = FALSE)
  h <- numeric(length = length(x))
  for(i in 1:length(h)){
    h[i] <- dweibull(x = x[i], shape = shape, scale = scale) / (1 - pweibull(q = x[i], shape = shape, scale = scale))
  }
  return(h)
}

# Generate values 
set.seed(12345)
weibull <- rweibull(n = 10000, shape = 1/2, scale = 1)
h <- wh(x = weibull, lambda = 1, kappa = 1/2)

# Plot the function 
pdf(file = "weibullhazards.pdf")

plot(x = c(1:length(h)), y = h, type = 'l', 
     main = "The Hazard Function of the Weibull", 
     xlab = "Time", ylab = "Hazard rate", ylim = c(0, 8), 
     bty = "n", col = "red", las = 1)

legend("topright", legend = c(TeX("$\\lambda = 1$, $\\kappa = 1/2$"), 
                           TeX("$\\lambda = 1$, $\\kappa = 1$"), 
                           TeX("$\\lambda = 1$, $\\kappa = 2$")), 
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

kappa <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(sigma)){
  weibull <- rweibull(n = 10000, shape = kappa[i], scale = 1)
  h <- wh(x = weibull, lambda = lambda[i], kappa = kappa[i])  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()

########### Gamma ########### (decreasing for shape < 1)

### Gamma density 

# Values to plot over 
x <- seq(from = 0, to = 20, length = 100000)

# Compute the densities 
d <- dgamma(x = x, shape = 0.9, rate = 1)

# Plot for different scale values
pdf(file = "gammadens.pdf")

plot(x, d, type = "l", lty = 1, xlab = "Values",
     ylab = "Density", main = "The Gamma Distribution", 
     bty = "n", ylim = c(0, 0.5), col = "red", las = 1)

legend("topright", legend = c(TeX("$\\alpha = 1$, $\\beta = 1$"), 
                              TeX("$\\alpha = 2$, $\\beta = 1$"), 
                              TeX("$\\alpha = 4$, $\\beta = 1$"), 
                              TeX("$\\alpha = 2$, $\\beta = 1/2$"), 
                              TeX("$\\alpha = 2$, $\\beta = 1/4$")),
       col = c("red", "blue", "dark green", "black", "orange"), lty = 1, cex = 0.8, bty = "n")

alpha <- c(2, 4, 2, 2)
beta <- c(1, 1, 1/2, 1/4)
col <- c("blue", "dark green", "black", "orange")

for (i in 1:length(alpha)){
  d <- dgamma(x = x, shape = alpha[i], rate = beta[i])
  lines(x, d, lwd = 1, col = col[i])
}

dev.off()

### Gamma hazard 
gammahaz <- function(x, shape, rate){
  x <- sort(x, decreasing = FALSE)
  h <- numeric(length = length(x))
  for(i in 1:length(h)){
    h[i] <- dgamma(x = x[i], shape, rate)/(1 - pgamma(q = x[i], shape, rate))
  }
  return(h)
}

## Compute hazards 
set.seed(12345)
x <- rgamma(n = 10000, shape = 2, rate = 1)
haz <- gammahaz(x, shape = 2, rate = 1)

# Plot the hazards 
pdf(file = "gammahazards.pdf")

plot(x = c(1:length(haz)), y = haz, type = 'l', 
     main = "The Hazard Function of the Gamma", 
     xlab = "Time", ylab = "Hazard rate", ylim = c(0, 5), 
     bty = "n", col = "red", las = 1)

legend("topright", legend = c(TeX("$\\alpha = 1/2$, $\\beta = 1$"), 
                              TeX("$\\alpha = 1$, $\\beta = 1$"), 
                              TeX("$\\alpha = 2$, $\\beta = 1$")), 
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

alpha <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(sigma)){
  gamma <- rgamma(n = 10000, shape = alpha[i], rate = 1)
  h <- gammahaz(x = gamma, shape = alpha[i], rate = 1)  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()

########## Fitting the network data to distributions ##########

### Import and visualize the results 
highnets <- readRDS(file = 'networks') # networks with high connectivity N = 1000
net = sort(highnets$time, decreasing = FALSE) 
hist(net, main = "Network hitting times")

## Hazard function of the network hitting times (monotonically decreasing)
data <- highnets 
data <- gbm1
data <- abm1
data <- std

surv_data <- data.frame(Time = data$time, Event = data$event, row.names = paste0("Sim", 1:length(data$time), ""))
surv_object <- Surv(time = data$time, event = data$event) 
surv_fit <- survfit(surv_object ~ 1, data = data, ctype = 1)
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data)

survival <- surv_fit$surv
surv_time <- surv_fit$time

plot(x = haz_fit$time, y = haz_fit$hazard, xlab = "Time", 
     ylab = "Hazard", ylim = c(min(haz_fit$hazard), max(haz_fit$hazard)),
     type = 'l', bty = 'n')#, main = 'Hazard function of GBM stopping times', col = 'blue')

## A function to check for how many points at the start the function increases 
trend <- function(x){
  trend <- logical(length = length(x-1))
  for(i in 2:length(x)){
    if(x[i-1] > x[i]){
      trend[i] <- TRUE
    } else {
      trend[i] <- FALSE
    }
  }
  return(trend)
}

trend <- trend(x = res)
length(which(trend == FALSE)) # how many increasing cases

### Test if the data come from Weibull (monotonically decreasing for k < 1)
hist(net)

### Test if the data come from lognormal (monotonically decreasing for μ > 1 and σ > 1)
hist(log(net)) # clearly not normal but may be due to the simulations. Gives an impression of normal. 
