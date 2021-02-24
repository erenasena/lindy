### Libraries 
library(statmod) # to work with the Inverse-Gaussian
library(VGAM) # to work with Levy 
library(fitdistrplus) # distribution fit 
library(parmsurvfit) # to fit right censored data 
library(parallel) # parallel computing for simulations 
library(latex2exp) # to write latex in plots 
library(survival) # survival analysis
library(survminer) # to plot the survival curves 
library(dynpred) # to compute conditional survival 
library(qualityTools) # for the Q-Q plots 
library(actuar) # for Pareto 

########################## Generating example hazard functions to display ########################## 
t <- c(1:100) # a time vector for all 

# The hazard function of Weibull 
wh <- function(x, lambda, k){
  x <- sort(x, decreasing = FALSE)
  h <- (k/lambda)*((x/lambda)^(k-1))
  return(h)
}

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
#pdf(file = "diff_hazards.pdf")
plot(x = t, y = w, type = "l", col = "red", lty = 1, ylim = c(0, 1), 
     xlab = "Time", ylab = "Hazard rate", main = "Hazard functions", las = 1, bty = "n")
lines(x = t, y = e, col = "blue", lty = 1, ylim = c(0, 1))
lines(x = t, y = p, col = "dark green", lty = 1, ylim = c(0, 1))
#dev.off()

########################## The Levy Distribution (hitting time of standard BM) ########################## 
n = 10000
#a = 0.1 # the distance between X0 and B 
#c = a^2

## Plotting the Levy density 

# Values to plot over 
x <- seq(from = 0.01, to = 3, length = n)

# Compute the densities 
d <- VGAM::dlevy(x = x, location = 0, scale = 1/2)

# Plot for different scale values
#pdf(file = "levydensities.pdf")

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

#dev.off()

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
#pdf(file = "levyhazards.pdf")

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

#dev.off()

########################## The Inverse Gaussian (IG) (distribution of ABM) ########################## 

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

########################## Lognormal ########################## decreasing for mu > 1, sigma > 1)

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

########################## Weibull ########################## decreasing for shape < 1)

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

legend("topright", legend = c(TeX("$\\lambda = 1$, $\\k = 1/2$"), 
                              TeX("$\\lambda = 1$, $\\k = 1$"), 
                              TeX("$\\lambda = 1$, $\\k = 3/2$"), 
                              TeX("$\\lambda = 1$, $\\k = 4$")), 
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
wh <- function(x, lambda, k){
  x <- sort(x, decreasing = FALSE)
  h <- (k/lambda)*((x/lambda)^(k-1))
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
h <- wh(x = weibull, lambda = 1, k = 1/2)

# Plot the function 
pdf(file = "weibullhazards.pdf")

plot(x = c(1:length(h)), y = h, type = 'l', 
     main = "The Hazard Function of the Weibull", 
     xlab = "Time", ylab = "Hazard rate", ylim = c(0, 8), 
     bty = "n", col = "red", las = 1)

legend("topright", legend = c(TeX("$\\lambda = 1$, $\\k = 1/2$"), 
                           TeX("$\\lambda = 1$, $\\k = 1$"), 
                           TeX("$\\lambda = 1$, $\\k = 2$")), 
       col = c("red", "blue", "dark green"), lty = 1, cex = 0.8, bty = "n")

k <- c(1, 2)
col <- c("blue", "dark green")

for (i in 1:length(k)){
  weibull <- rweibull(n = 10000, shape = k[i], scale = 1)
  h <- wh(x = weibull, lambda = 1, k = k[i])  
  lines(x = c(1:length(h)), y = h, lwd = 1, col = col[i])
}

dev.off()

########################## Gamma ########################## decreasing for shape < 1)

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

########################## Fitting the network data to distributions ########################## 

### Import and visualize the results 
highnets <- readRDS(file = 'networks') # networks with high connectivity N = 1000
highnets = sort(highnets$time, decreasing = FALSE) 
lownets <- readRDS(file = 'lownets')
lownets <- sort(lownets$time, decreasing = FALSE)

### Distribution fit with MLE  

## Prepare the data for right-censored estimation
cens <- function(data){ # network lifetimes
  cens <- numeric(length = length(data))
  for(i in 1:length(cens)){
    if(data[i] >= max(data)){
      cens[i] <- NA
    }else if(data[i] < max(data)){
      cens[i] <- data[i]
    }
  }
  return(cens)
}

nets <- lownets
cens <- cens(data = nets)
length(which(nets >= max(nets))) == sum(is.na(cens)) # testing if it worked 
data <- data.frame("left" = nets, "right" = cens)

# Fit distributions
fit.weibull <- fitdistcens(censdata = data, distr = "weibull")
fit.lognormal <- fitdistcens(censdata = data, distr = "lnorm")
fit.gamma <- fitdistcens(censdata = data, distr = "gamma", lower = c(0, 0))
fit.exp <- fitdistcens(censdata = data, distr = "exp")

## Uncertainty estimation 
detectCores(all.tests = FALSE, logical = TRUE) # detect n of cores for parallel processing

# Uncertainty in the parameters

# Weibull 
weibull_params <- bootdistcens(f = fit.weibull, niter = 1000, silent = FALSE, 
                               parallel = "multicore", ncpus = 8)

weibull_params$fitpart
weibull_params$CI

# Lognormal 
lnorm_params <- bootdistcens(f = fit.lognormal, niter = 1000, silent = FALSE, 
                             parallel = "multicore", ncpus = 8)

lnorm_params$fitpart
lnorm_params$CI

# Gamma 
gamma_params <- bootdistcens(f = fit.gamma, niter = 1000, silent = FALSE, 
                             parallel = "multicore", ncpus = 8)

gamma_params$fitpart
gamma_params$CI

# Exponential 
exp_params <- bootdistcens(f = fit.exp, niter = 1000, silent = FALSE, 
                             parallel = "multicore", ncpus = 8)

exp_params$fitpart
exp_params$CI

# Uncertainty in the CDF 

# Weibull 
CIcdfplot(weibull_params, CI.level= 0.95, CI.output = "probability", 
          main = "Empirical and theoretical CDF for Weibull", 
          xlab = "Censored lifetimes", ylab = TeX('$\\P(X \\leq x)$'), 
          las = 1, bty = "n")

legend("bottomright", legend = c("Empirical", "CI", "Theoretical"),
       col = c("black", "red", "red"), lty = 1:2, cex = 1, bty = "n")


# Lognormal 
CIcdfplot(lnorm_params, CI.level= 0.95, CI.output = "probability", 
          main = "Empirical and theoretical CDF for Lognormal", 
          xlab = "Censored lifetimes", ylab = TeX('$\\P(X \\leq x)$'), 
          las = 1, bty = "n")

legend("bottomright", legend = c("Empirical", "CI", "Theoretical"),
       col = c("black", "red", "red"), lty = 1:2, cex = 1, bty = "n")


# Gamma 
CIcdfplot(gamma_params, CI.level= 0.95, CI.output = "probability", 
          main = "Empirical and theoretical CDF for Gamma", 
          xlab = "Censored lifetimes", ylab = TeX('$\\P(X \\leq x)$'), 
          las = 1, bty = "n")

legend("bottomright", legend = c("Empirical", "CI", "Theoretical"),
       col = c("black", "red", "red"), lty = 1:2, cex = 1, bty = "n")

# Exponential 
CIcdfplot(exp_params, CI.level= 0.95, CI.output = "probability", 
          main = "Empirical and theoretical CDF for Exponential", 
          xlab = "Censored lifetimes", ylab = TeX('$\\P(X \\leq x)$'), 
          las = 1, bty = "n")

legend("bottomright", legend = c("Empirical", "CI", "Theoretical"),
       col = c("black", "red", "red"), lty = 1:2, cex = 1, bty = "n")

## Absolute fit
censor <- numeric(length = length(nets))
for(i in 1:length(censor)){
  if(nets[i] >= max(nets)){
    censor[i] <- 0
  } else if(nets[i] < max(nets)){
    censor[i] <- 1
  }
}

data <- data.frame("time" = nets, "censor" = censor)

# Weibull 
AD_weibull <- compute_AD(data, dist = "weibull", time = "time", censor = "censor")
print(AD_weibull)

# Lognormal 
AD_lognormal <- compute_AD(data, dist = "lnorm", time = "time", censor = "censor")
print(AD_lognormal)

# Exponential
AD_exp <- compute_AD(data, dist = "exp", time = "time", censor = "censor")
print(AD_exp)

# Gamma 
AD_gamma <- compute_AD(data, dist = "gamma", time = "time", censor = "censor")
print(AD_gamma)

## Relative fit (good stuff)
summary(fit.weibull)
summary(fit.lognormal)
summary(fit.gamma)
summary(fit.exp)

(fit.lognormal$aic < fit.weibull$aic) & (fit.weibull$aic < fit.gamma$aic)
(fit.lognormal$bic < fit.weibull$bic) & (fit.weibull$bic < fit.gamma$bic)

fit.lognormal$aic - fit.weibull$aic # how big is the difference between lognormal and Weibull? 

# By plotting the ecdf
l <- list(fit.weibull, fit.gamma)

pdf(file = "weakecdfvscdf.pdf")

cdfcompcens(ft = l, xlogscale = FALSE, ylogscale = FALSE, 
            main = "Comparing ECDF of the networks with theoretical CDFs",
            xlab = "Censored lifetimes", ylab = TeX('$\\P(X \\leq x)$'), 
            las = 1, bty = "n") 

legend("bottomright", legend = c("Weibull", "Gamma"),
       col = c("red", "green"), lty = 1:2, cex = 1, bty = "n")

dev.off()

# By plotting the Q-Q plot 
pdf(file = "weakqqcomparison.pdf")

qqcompcens(ft = l, xlogscale = FALSE, ylogscale = FALSE, 
           main = "Q-Q plot of Weibull and Gamma distributions", 
           xlab = "Theoretical quantiles", ylab = "Empirical quantiles", 
           las = 1, bty = "n") 

legend("bottomright", legend = c("Weibull", "Gamma"),
       col = c("red", "green"), lty = 1:2, cex = 1, bty = "n")

dev.off()
?qqcompcens

########################## Original stuff ########################## 
### Load the data 
highnets <- readRDS(file = 'networks') # networks with high connectivity 
lownets <- readRDS(file = 'lownets') # networks with low connectivity 

### Visualizing the distributions

# Some distributions
set.seed(12345)
g <- rgamma(n = 100000, shape = 0.5, rate = 1)
w <- rweibull(n = 100000, shape = 0.9, scale = 50)
e <- rexp(n = 100000, 1)
p <- rpareto(n = 100000, scale = 1, shape = 3)

# Histogram of the hitting times
pdf(file = "stronghist.pdf")
hist(highnets$time, breaks = 25, xlim = c(0, 10000), main = "Frequency distribution of 
     strongly connected networks", xlab = "Time", col = "lightblue", border = "darkblue", prob = F)
dev.off()

pdf(file = "weakhist.pdf")
hist(lownets$time, breaks = 25, xlim = c(0, 70), main = "Frequency distribution of 
     weakly connected networks", xlab = "Time", col = "lightblue", border = "darkblue", prob = F)
dev.off()

# Q-Q Plot to check exponentiality (if linear, thin tails, if concave, there may be heavy tailedness)
hittings <- sort(net) # sort the data 
p <- ppoints(hittings, length(hittings)) # get the probabilities of the data
s <- quantile(x = hittings, p = p) # sample quantiles
q <- qexp(p = p) # exponential quantiles 
qqplot(x = s, y = q, xlab = "Sample quantiles", ylab = "Theoretical quantiles", 
       main = "Exponential QQ Plot")

png(file = "qqhigh.png")
qualityTools::qqPlot(x = net, y = "exponential", confbounds = TRUE, 
                     bty = "n", main = "Exponential Q-Q plot of highly connected networks")
dev.off()

png(file = "qqlow.png")
qualityTools::qqPlot(x = lownets, y = "exponential", confbounds = TRUE, 
                     bty = "n", main = "Exponential Q-Q plot of weakly connected networks")
dev.off()

# The Zipf plot function of Cirillo 
zipfplot <- function (data, type = "plot", title = TRUE){
  # type should be equal to ’points’ if you want to add the
  # Zipf Plot to an existing graph
  # With other strings or no string a new graph is created.
  # If title is set to be F, the title of the plot is not given.
  # This can be useful when embedding the Zipf plot into other
  # plots.
  data <- sort(as.numeric(data)) #sorting data 
  y <- 1 - ppoints(data) #computing 1-F(x)
  if (type == "points") {
    points(data, y, xlog = T, ylog = T, xlab = "x on log scale",
           ylab = "1-F(x) on log scale")
  }
  plot(data, y, log = "xy", xlab = "x on log scale", ylab = "1-F(x) on log scale", 
       main = "Zipf Plot of highly connected networks", bty = "n")
}

#png(filename = "zipfcens2.png")
zipfplot(highnets$time)
#dev.off()

#png(filename = "lowzipf.png")
zipfplot(lownets$time)
#dev.off()

#png(filename = "zipfweibull")
zipfplot(data = w)
#dev.off()

#png(filename = "zipfpareto")
zipfplot(data = p)
#dev.off()

########################## Survival analysis ########################## 
?bshazard

## Create survival objects

# Highly connected networks 
surv_data <- data.frame(Time = highnets$time, Event = highnets$event, 
                        row.names = paste0("Sim", 1:length(highnets$time), ""))
surv_object <- Surv(time = highnets$time, event = highnets$event) 
surv_fit <- survfit(surv_object ~ 1, data = highnets, ctype = 1)

# Weakly connected networks 
surv_data2 <- data.frame(Time = lownets$time, Event = lownets$event, 
                         row.names = paste0("Sim", 1:length(lownets$time), ""))
surv_object2 <- Surv(time = lownets$time, event = lownets$event) 
surv_fit2 <- survfit(surv_object2 ~ 1, data = lownets, ctype = 1)

# Highly connected networks 
haz_fit <- bshazard::bshazard(surv_object ~ 1, data = surv_data, alpha = 0.05)
haz_fit$phi # overdispersion parameter
summary(haz_fit)
?bshazard

# Weakly connected networks 
haz_fit2 <- bshazard::bshazard(surv_object2 ~ 1, data = surv_data2, alpha = 0.05)
haz_fit2$phi

## Plot the results

# Highly connected networks
pdf(file = "stronghazard.pdf")

plot(x = haz_fit$time, y = haz_fit$hazard, xlab = "Time", 
     ylab = "Hazard", ylim = c(min(haz_fit$hazard), max(haz_fit$hazard)),
     type = 'l', bty = 'n', main = 'Hazard function of strongly connected networks')
points(haz_fit$raw.data$time, haz_fit$raw.data$raw.hazard, cex = .1, lwd = 3, col = 1)
lines(haz_fit$time, haz_fit$hazard, lty = 1, lwd = 2)
lines(haz_fit$time, haz_fit$lower.ci, lty = 3, lwd = 1, col = "red")
lines(haz_fit$time, haz_fit$upper.ci, lty = 3, lwd = 1, col = "red")

dev.off()

# Weakly connected networks
pdf(file = "weakhazard.pdf")

plot(x = haz_fit2$time, y = haz_fit2$hazard, xlab = "Time", 
     ylab = "Hazard", ylim = c(min(haz_fit2$hazard), max(haz_fit2$hazard)),
     type = 'l', bty = 'n', main = 'Hazard function of weakly connected networks')
points(haz_fit2$raw.data$time, haz_fit2$raw.data$raw.hazard, cex = .1, lwd = 3, col = 1)
lines(haz_fit2$time, haz_fit2$hazard, lty = 1, lwd = 2)
lines(haz_fit2$time, haz_fit2$lower.ci, lty = 3, lwd = 1, col = "red")
lines(haz_fit2$time, haz_fit2$upper.ci, lty = 3, lwd = 1, col = "red")

dev.off()

########################## Conditional Survival ########################## 
#t = c(1:10000)

### Pareto 
#cspareto <- function(t, s, alpha){
  #cs <- numeric(length = length(t))
  #for(i in 1:length(cs)){
    #cs[i] <- (t[i]/(t[i]+s))^alpha
  #}
  #return(cs)
#}

#cs <- cspareto(t = t, s = 1000, alpha = 2)
#plot(x = t, y = cs, type = "l")

### From a sample 

## Strongly connected networks
fit <- dynpred::Fwindow(object = surv_fit, width = 10, variance = TRUE, conf.level = 0.95)
condeath <- fit$Fw
consurv <- 1 - condeath

#pdf(file = "strongcond.pdf")

plot(x = fit$time, y = consurv, type = 'l', col = 'green', xlab = 'Time', 
ylab = 'Conditional Probability', main = 'Conditional survival and death probabilities of 
strongly connected networks',
ylim = c(0.95, 1), xlim = c(0, 10000), bty = 'n')
lines(x = fit$time, y = fit$low)
lines(x = fit$time, y = fit$up)

lines(x = fit$time, y = condeath, col = 'red') # change the limit of the y-axis to c(0, 1) to see this
lines(x = fit$time, y = 1-fit$low)
lines(x = fit$time, y = 1-fit$up)

legend(x = 'right', bty = 'n', lty = c(1, 1), col = c("green", "red"), 
       legend = c("Conditional survival", "Conditional death"))

#dev.off()

## Weakly connected networks
fit2 <- dynpred::Fwindow(object = surv_fit2, width = 10, variance = TRUE, conf.level = 0.95)
condeath2 <- fit2$Fw
consurv2 <- 1 - condeath2

pdf(file = "weakconds.pdf")

plot(x = fit2$time, y = consurv2, type = 'l', col = 'green', xlab = 'Time', 
     ylab = 'Conditional Probability', main = 'Conditional survival and death probabilities of 
weakly connected networks', bty = 'n')
lines(x = fit2$time, y = fit2$low)
lines(x = fit2$time, y = fit2$up)

#lines(x = fit2$time, y = condeath2, col = 'red') # change the limit of the y-axis to c(0, 1) to see this
#lines(x = fit2$time, y = 1-fit2$low)
#lines(x = fit2$time, y = 1-fit2$up)

legend(x = 'right', bty = 'n', lty = c(1, 1), col = "green", legend = "Conditional survival")

dev.off()

#cond <- data.frame(con_time, con_surv, con_death)
#colnames(cond) <- c('Time', 'Conditional Survival', 'Conditional Death')