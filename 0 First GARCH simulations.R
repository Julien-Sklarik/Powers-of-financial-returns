#package needed for GARCH model
#install.packages("tseries")
require("tseries")
library(tidyverse)
library(ggthemes)
library(forecast)
library(tseries)
library(gridExtra)
library(rugarch)

set.seed(1)

#the coefficients for our GARCH(1,1) model
w <- 1 #R nota: a0
# I'm choosing alpha and beta as in the book, but it's also because an alpha small and a beta large are typical of those of stock returns.
alpha <- 0.2 #R nota: a1
beta <- 0.7 #R nota: b1
# I want to have alpha + beta < 1, because this way it will exist a unique stationary solution. 
# The variance will be stationary.
# The closer the sum alpha + beta is to 1, the slower is the decay in the autocorrelation function of squared residuals.

#For eta, I take the following gaussian strong white noise term values (which is a white noise of variance = 1)
eta <- rnorm(2000)

#the actual x(t) time series
x <-rep(0,2000)
#volatility squared values
sigma2 <- rep(0,2000)

#GARCH(1,1) model simulation
for(t in 2:2000){
  sigma2[t] <- w+alpha*(x[t-1]^2)+beta*sigma2[t-1]
  x[t] <- eta[t]*sqrt(sigma2[t])
}

plot(x)
#autocorrelation plots
acf(x)
pacf(x)
acf(x*x) # Here we could check the 6th stylised fact of Volatility clustering, if we were using real data.
# My ref for the stylised facts : http://rama.cont.perso.math.cnrs.fr/pdf/empirical.pdf

#use the GARCH function
x.garch <- garch(x,trace=FALSE)
#show the confidence intervals for the parameters
confint(x.garch)
# 2.5 %     97.5 %
# a0  8.1720073 12.8800103 ## 1 is far from this interval
# a1  0.1201616  0.2223641
# b1 -0.1572084  0.1572084

# Here I try to apply automatically the Box-Jenkins methodology to see which ARIMA is proposed.
# Box-Jenkins approach applies ARIMA models to find the best fit of a time series model that represent the stochastic process which generate time series. This method uses a three stage modelling approach: a) identification, b) estimation, c) diagnostic checking.

model.arima = forecast::auto.arima(x , max.order = c(3 , 0 ,3) , stationary = TRUE , trace = T , ic = 'aicc')
# Best model: ARIMA(0,0,0) with zero mean













