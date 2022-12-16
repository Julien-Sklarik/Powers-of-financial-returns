rm(list = ls())
.rs.restartR()
#package needed for GARCH model
#install.packages("tseries")
#install.packages("quantmod") 
library(quantmod)  
require("tseries")
set.seed(1)

#the coefficients for our GARCH(1,1) model
w <- 1 #R nota: a0
# I'm choosing alpha and beta as in the book, but it's also because an alpha small and a beta large are typical of those of stock returns.
alpha <- 0.2 #R nota: a1
beta <- 0.7 #R nota: b1
# I want to have alpha + beta < 1, because this way it will exist a unique stationary solution. 
# The variance will be stationary.
# The closer the sum alpha + beta is to 1, the slower is the decay in the autocorrelation function of squared residuals.
T <- 5000 # The number of obs
#For eta, I take the following gaussian strong white noise term values (which is a white noise of variance = 1)
eta <- rnorm(T)
#the actual x(t) time series
x <-rep(0,T)

#volatility squared values
sigma2 <- rep(0,T)

#GARCH(1,1) model simulation
for(t in 2:T){
  sigma2[t] <- w+alpha*(x[t-1]^2)+beta*sigma2[t-1]
  x[t] <- eta[t]*sqrt(sigma2[t])
}
plot(x, type = "l")

# # Here is just a code to share my garch11 with the other, so that they can guess the parameters with their functions
# process <- x
# save(process, file= "/Users/juliensklarik/Documents/Stat app/Data/process.RData") # alpha = 0.14 and beta = 0.65
# load("/Users/juliensklarik/Documents/Stat app/Data/process.RData")
# write.csv(process, "/Users/juliensklarik/Documents/Stat app/Data/process.csv", row.names=FALSE)
# # Reading Georgii's data:
# variables <- read.csv(file = '/Users/juliensklarik/Documents/Stat app/Data/garch_9000_values_georgii.csv')
# x <- variables$X0


garch11Fit = function(x) {
  require(numDeriv)
  ### Initialization of Model Parameters and Bounds:
  # we initialize the set of model parameters {θ}, params, and the corresponding upper and lower bounds. 
  # We use bounds lowerBounds and upperBounds which are wide enough to serve almost every economic and financial GARCH(1,1) model, 
  # and define model parameters which take typical values.
  Mean = mean(x)
  Var = var(x)
  S = 1e-06
  params = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8) #  We initialize omega by the variance of the series adjusted by the persistence omega = Var(x) ∗ (1 − a − b)
  lowerBounds = c(mu = -10 * abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10 * abs(Mean), omega = 100 * Var, alpha = 1 - S, beta = 1 - S)
  ### Set Conditional Distribution Function:
  # For the conditional distribution we use the Normal distribution dnorm().
  garchDist = function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  ### Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (x - mu)
    Mean = mean(z^2) # The variance of x
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2) # The vector c(Mean, z[-length(x)]^2) is z^2 where we suppress the last value and add Mean at the beginning. 
    # We do so because the last value won't be use. Each sigma^2 is calculated with t-1 values. sigma^2(t=1) is with epsilon_0^2= mean^2=sigma_0^2 (I saw that it was usually the choice of authors, but I don't know why. It looks like an arbitrary choice, to some extent.)
    # above, we are preparing the computation of sigma^2= omega+alpha*x_(t-1)^2+beta*sigma_(t-1)^2, the last part is obtained below
    h = filter(e, beta, "r", init = Mean) # r for recurssive: an autoregression is used.
    # h correspond to the series of sigma^2 values!!
    hh = sqrt(abs(h)) # as h are empirical values, we use abs to correct the potential negative variance estimations
    -sum(log(garchDist(z, hh)))  #here we have exactly the negative log-likelihood for the GARCH(1,1) model
  }
  # We could add the following line:
  # print(garchLLH(params)) 
  
  ### Estimate Parameters and Compute Numerically the Hessian:
  # For the GARCH(1,1) model optimization of the log-likelihood function we use the constrained solver nlminb(). 
  # To compute standard errors and t-values we evaluate the Hessian matrix numerically.
  fit = nlminb(start = params, objective = garchLLH, lower = lowerBounds, upper = upperBounds) # To avoid using those functions, I could code my own maximiser function. But it's not the objective here.
  Hessian <- numDeriv::hessian(func = garchLLH, x = fit$par) # Hessian is just a matrix 4x4, bc 4 is the number of parameters
  
  ### Create and Print Summary Report:
  # The results for the estimated parameters together with standard errors and t-values are summarized and printed.
  se.coef = sqrt(diag(solve(Hessian))) # solve gives the matrix 4x4 X such that H %*% X == I
  # diag gives the diagonal vector and we are obtaining the s.e. bc the variance is the inverse of Fisher Information, which is -E(H(theta)) by definition.
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval)))) # We are construction our dataframe
  dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)")) # We rename our columns and rows
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
}

# Fitting a GARCH(1,1) model with normal errors to the TS.
garch11_model <- garch11Fit(x)
# Both a1 and β1 are significantly different from zero, therefore it is reasonable to assume time-varying volatility of the residuals.



## Version of the normal garch fit function which returns only the three coeff of interest: w, alpha and beta

garch11Fit.coef = function(x) {
  require(numDeriv)
  ### Initialization of Model Parameters and Bounds:
  # we initialize the set of model parameters {θ}, params, and the corresponding upper and lower bounds. 
  # We use bounds lowerBounds and upperBounds which are wide enough to serve almost every economic and financial GARCH(1,1) model, 
  # and define model parameters which take typical values.
  Mean = mean(x)
  Var = var(x)
  S = 1e-06
  params = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8) #  We initialize omega by the variance of the series adjusted by the persistence omega = Var(x) ∗ (1 − a − b)
  lowerBounds = c(mu = -10 * abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10 * abs(Mean), omega = 100 * Var, alpha = 1 - S, beta = 1 - S)
  ### Set Conditional Distribution Function:
  # For the conditional distribution we use the Normal distribution dnorm().
  garchDist = function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  ### Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (x - mu)
    Mean = mean(z^2) # The variance of x
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2) # The vector c(Mean, z[-length(x)]^2) is z^2 where we suppress the last value and add Mean at the beginning. 
    # We do so because the last value won't be use. Each sigma^2 is calculated with t-1 values. sigma^2(t=1) is with epsilon_0^2= mean^2=sigma_0^2 (I saw that it was the choice of authors, but I don't know why. It looks like an arbitrary choice, to some extent.)
    # above, we are preparing the computation of sigma^2= omega+alpha*x_(t-1)^2+beta*sigma_(t-1)^2, the last part is obtained below
    h = filter(e, beta, "r", init = Mean) # r for recurssive: an autoregression is used.
    # h correspond to the series of sigma^2 values!!
    hh = sqrt(abs(h)) # as h are empirical values, we use abs to correct the potential negative variance estimations
    -sum(log(garchDist(z, hh)))  #here we have exactly the negative log-likelihood for the GARCH(1,1) model
  }
  # We could add the following line:
  # print(garchLLH(params)) 
  
  ### Estimate Parameters and Compute Numerically the Hessian:
  # For the GARCH(1,1) model optimization of the log-likelihood function we use the constrained solver nlminb(). 
  # To compute standard errors and t-values we evaluate the Hessian matrix numerically.
  fit = nlminb(start = params, objective = garchLLH, lower = lowerBounds, upper = upperBounds) # To avoid using those functions, I could code my own maximiser function. But it's not the objective here.
  fit$par[-1]
}

garch11_model <- rbind(garch11Fit.coef(x),garch11Fit.coef(x))
garch11_model[,1]

N <- 1000 # The number of simulations of size T (or n)
T <-500 # The number of observations (We could call it n.)
w <- 1 #R nota: a0
alpha <- 0.2 #R nota: a1
beta <- 0.7 #R nota: b1

validation.garch11Fit.coef = function(N,T) {
  coef <- c()
  x <-rep(0,T)
  sigma2 <- rep(0,T)
  for (k in 1:(N-1)) {
    #GARCH(1,1) model simulation
    eta <- rnorm(T)
    for(t in 2:T){
      sigma2[t] <- w+alpha*(x[t-1]^2)+beta*sigma2[t-1]
      x[t] <- eta[t]*sqrt(sigma2[t])
    }
    coef <- rbind(garch11Fit.coef(x)-c(w,alpha,beta), coef)
  }
  boxplot(coef)
  coef
}

validation.garch11Fit.coef(N,T)



## Variance of the coeffs

coeffs.garch11Fit.coef = function(N,T) {
  coef <- c()
  x <-rep(0,T)
  sigma2 <- rep(0,T)
  for (k in 1:(N-1)) {
    #GARCH(1,1) model simulation
    eta <- rnorm(T)
    for(t in 2:T){
      sigma2[t] <- w+alpha*(x[t-1]^2)+beta*sigma2[t-1]
      x[t] <- eta[t]*sqrt(sigma2[t])
    }
    coef <- rbind(garch11Fit.coef(x), coef)
  }
  coef
}
coeffs <- coeffs.garch11Fit.coef(N,T)
var(coeffs)



# Computation of J^-1 and comparison with the hessian computed in garch11fit:

j = function(x){
  Mean = mean(x)
  Var = var(x)
  S = 1e-06
  n = length(x)
  parm = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8)
  mu = parm[1]
  omega = parm[2]
  alpha = parm[3]
  beta = parm[4]
  z = (x - mu)
  Mean = mean(z^2) # The variance of x
  # Use Filter Representation:
  e = omega + alpha * c(Mean, z[-length(x)]^2) # The vector c(Mean, z[-length(x)]^2) is z^2 where we suppress the last value and add Mean at the beginning. 
  # We do so because the last value won't be use. Each sigma^2 is calculated with t-1 values. sigma^2(t=1) is with epsilon_0^2= mean^2=sigma_0^2 (I saw that it was the choice of authors, but I don't know why. It looks like an arbitrary choice, to some extent.)
  # above, we are preparing the computation of sigma^2= omega+alpha*x_(t-1)^2+beta*sigma_(t-1)^2, the last part is obtained below
  h = filter(e, beta, "r", init = Mean) # r for recurssive: an autoregression is used.
  # h correspond to the series of sigma^2 values!!
  hh = sqrt(abs(h))
  
  mat1.data <- rep(0,9)
  j.sum <- matrix(mat1.data,nrow=3,ncol=3,byrow=TRUE)
  s2w <- 1/(1-beta)
  for (t in 2:n){
    s2a <- 0
    s2b <- 0
    for (k in 1:(t-1)){
      s2a <- s2a + (beta^(k+1))*z[t-k]^2
      s2b <- s2b + (beta^(k+1))*h[t-k]
    }
    s2t <- matrix(c(s2w, s2a, s2b),nrow=3,ncol=1,byrow=TRUE)
    s2tp <- matrix(c(s2w, s2a, s2b),nrow=1,ncol=3,byrow=TRUE)
    j.sum <- j.sum + (1/h[t])*(s2t %*% s2tp) #here I must multiply the two vector to make a matrix
  }
  print(solve(j.sum/(n-1)))
  
  garchDist = function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  lowerBounds = c(mu = -10 * abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10 * abs(Mean), omega = 100 * Var, alpha = 1 - S, beta = 1 - S)
  garchLLH = function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (x - mu)
    Mean = mean(z^2) # The variance of x
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2) # The vector c(Mean, z[-length(x)]^2) is z^2 where we suppress the last value and add Mean at the beginning. 
    # We do so because the last value won't be use. Each sigma^2 is calculated with t-1 values. sigma^2(t=1) is with epsilon_0^2= mean^2=sigma_0^2 (I saw that it was the choice of authors, but I don't know why. It looks like an arbitrary choice, to some extent.)
    # above, we are preparing the computation of sigma^2= omega+alpha*x_(t-1)^2+beta*sigma_(t-1)^2, the last part is obtained below
    h = filter(e, beta, "r", init = Mean) # r for recurssive: an autoregression is used.
    # h correspond to the series of sigma^2 values!!
    hh = sqrt(abs(h)) # as h are empirical values, we use abs to correct the potential negative variance estimations
    -sum(log(garchDist(z, hh)))  #here we have exactly the negative log-likelihood for the GARCH(1,1) model
  }
  # We could add the following line:
  # print(garchLLH(params)) 
  
  ### Estimate Parameters and Compute Numerically the Hessian:
  # For the GARCH(1,1) model optimization of the log-likelihood function we use the constrained solver nlminb(). 
  # To compute standard errors and t-values we evaluate the Hessian matrix numerically.
  fit = nlminb(start = parm, objective = garchLLH, lower = lowerBounds, upper = upperBounds) # To avoid using those functions, I could code my own maximiser function. But it's not the objective here.
  Hessian <- numDeriv::hessian(func = garchLLH, x = fit$par) # Hessian is just a matrix 4x4, bc 4 is the number of parameters
  Hessian
}
j(x)


























### Trials and failures:
Mean = mean(x)
Var = var(x)
S = 1e-06
parm = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8)
mu = parm[1]
omega = parm[2]
alpha = parm[3]
beta = parm[4]
z = (x - mu)
Mean = mean(z^2)
# Use Filter Representation:
e = omega + alpha * c(Mean, z[-length(x)]^2) # 
h = filter(e, beta, "r", init = Mean) # r for recurssive: an autoregression is used.
hh = sqrt(abs(h))
garchDist = function(z, hh) {
  dnorm(x = z/hh)/hh
}
-sum(log(garchDist(z, hh)))  

filter(c(1,2,3,4), 2, "r", init = 1)






































