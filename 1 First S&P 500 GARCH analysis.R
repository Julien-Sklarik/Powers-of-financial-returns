# The package quantmod stands for Quantitative Financial Modeling Framework. 
# This allows to web-scrape financial data.
#install.packages("quantmod") 
library(quantmod)            

getSymbols(Symbols = "^GSPC", 
           src = "yahoo",     
           from="2003-01-01",
           to="2021-12-31",
           verbose=F)

chartSeries(GSPC, theme="white")

data(GSPC$Close)
plot.ts(nottem)


# 2nd trial to have a better format, with the function get.hist.quote of the package tseries.
library(rugarch)
library(tseries)
gspc = get.hist.quote(instrument = "IBM", quote = "Adj", provider = c("yahoo"), 
                     method = NULL, start = "2003-01-01", end = "2021-12-31", compression = "d", 
                     retclass = c("zoo"), quiet = FALSE, drop = FALSE)

#Transform to log-returns
gspc_prices <- as.data.frame(gspc)
N <- length(gspc_prices[, 1])
gspc_prices.returns <- 100*(log(gspc_prices[2:N, ])-log(gspc_prices[1:(N-1), ]))

# ACF and PACF of the Squared residuals
TSA::acf(I(residuals(lm(gspc_prices.returns ~ 1))^2), main = "Squared residuals")

pacf(I(residuals(lm(gspc_prices.returns ~ 1))^2), main = "Squared residuals")

# Fitting a GARCH(1,1) model with normal errors to the TS.
garch11_model <- garch11Fit(gspc_prices.returns)
# Both a1 and β1 are significantly different from zero, therefore it is reasonable to assume time-varying volatility of the residuals.

# How does one proceed with the estimation of a GARCH model? 
# Maximum likelihood is the standard option, but the MLE must be found numerically. 
# This function from a preprint by Würtz, Chalabi and Luskan, shows how to construct the likelihood for a simple GARCH(1,1) model.
# Here is the reference :
# http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf 
# From Wurtz, Chalabi, Luskan (JSS)
garch11Fit = function(x) {
  require(numDeriv)
  # Step 1: Initialize Model Parameters and Bounds:
  Mean = mean(x)
  Var = var(x)
  S = 1e-06
  params = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8)
  lowerBounds = c(mu = -10 * abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10 * abs(Mean), omega = 100 * Var, alpha = 1 - S, beta = 1 - 
                    S)
  # Step 2: Set Conditional Distribution Function:
  garchDist = function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  # Step 3: Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (x - mu)
    Mean = mean(z^2)
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2)
    h = filter(e, beta, "r", init = Mean)
    hh = sqrt(abs(h))
    -sum(log(garchDist(z, hh)))  #llh
    
  }
  # print(garchLLH(params)) Step 4: Estimate Parameters and Compute
  # Numerically Hessian:
  fit = nlminb(start = params, objective = garchLLH, lower = lowerBounds, 
               upper = upperBounds)
  Hessian <- numDeriv::hessian(func = garchLLH, x = fit$par)
  # Step 5: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", 
                                          "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
}


garch11_model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                            mean.model = list(include.mean = TRUE), 
                            distribution.model = "norm")
# Model fitting
garch11_fit <- ugarchfit(spec = garch11_model, data = gspc_lret, 
                         solver.control = list(solver = 9), solver ='hybrid')
qqnorm(scale(residuals(garch11_fit)), pty = "s")
abline(a = 0, b = 1)

sres <- residuals(garch11_fit)/sigma(garch11_fit)
# residuals(garch11_fit, standardize = TRUE) sres^2 is equivalent to
# residuals(garch11_fit)^2/garch11_fit@fit$var
TSA::acf(sres^2, main = "Squared standardized residuals of GARCH(1,1) model")
pacf(sres^2, main = "Squared standardized residuals of GARCH(1,1) model")

# Parameters of the GARCH(1,1)
alpha1 <- garch11_fit@model$pars["alpha1", 1]
beta1 <- garch11_fit@model$pars["beta1", 1]

print(garch11_fit)

# kurtosis_garch11 <- function(alpha1, beta1) {
#   3 * (1 + alpha1 + beta1) * (1 - alpha1 - beta1)/(1 - beta1^2 - 2 * alpha1 * 
#                                                      beta1 - 3 * alpha1^2)
# }
# kurtosis_garch11(garch11_model[3, 1], garch11_model[4, 1])

# We see that we are close to having an IGARCH as alpha + beta is close to 1.  
# The fact that beta is close to one is typical of heavy-tailed financial returns. 
# The kurtosis is provided if the fourth moment exist and 3a^2 + 2ab + b^2<1 
# Thus: A shock at one time becomes permanent, the conditional variance is not stationary.


# Just to test, I will try a GARCH(2,1):
library(rugarch)
# Higher order needed (?), also more heavy-tailed innovations
garch21_model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 
                                                                                   1)), mean.model = list(include.mean = TRUE), distribution.model = "std")
garch21_fit <- ugarchfit(spec = garch21_model, data = ibm_lret, solver = "hybrid", 
                         solver.control = list(solver = 9))
# Is the conditional variance model stationary?
garch21_fit@model$pars["alpha1", 1] + garch21_fit@model$pars["alpha2", 1] + 
  garch21_fit@model$pars["beta1", 1]
# Yes
plot(garch21_fit, which = 9)
plot(garch21_fit, which = 11)
# TSA::acf(c((residuals(garch21_fit)/sigma(garch21_fit))^2), main = 'Scaled
# residuals') pacf(c((residuals(garch21_fit)/sigma(garch21_fit))^2), main =
# 'Scaled residuals')
gjrgarch_model <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 
                                                                                      1)), mean.model = list(include.mean = TRUE), distribution.model = "std")
gjrgarch_fit <- ugarchfit(spec = gjrgarch_model, data = ibm_lret, solver = "hybrid", 
                          solver.control = list(solver = 9))
# Well, I don't think it improve anything, but it was worth a try ^^'






















######################### Failed trials...

#The following commands will compute GARCH(m,s). Keep in mind that it may not converge for certain combinations of m and s.

Apple_garch <-  ugarchspec(variance.model = list(model="sGARCH",         #Other options are egarch, fgarch, etc.
                                                 garchOrder=c(1,2)), # You can modify the order GARCH(m,s) here
                           mean.model = list(armaOrder=c(3,2)), #Specify your ARMA model implying your model should be stationary.
                           distribution.model = "norm")         #Other distribution are "std" for t-distribution, and "ged" for General Error Distribution
Apple_Garch2 <- ugarchfit(spec=Apple_garch, 
                          data=gspc_prices.returns)

garch <- Apple_Garch2@sigma.t   #this is your volatility series


install.packages("ggplot2")
library(ggplot2)
garch_plot <- ggplot2()+
  geom_line(aes(x=as.numeric(row.names(gspc_prices.returns)),
                y=Apple,
                color="log returns"),
            data=gspc_prices.returns)+
  geom_line(aes(x=as.numeric(row.names(gspc_prices.returns)),
                y=Apple_Garch2,
                color="Volatility"))+
  geom_line(aes(x=as.numeric(row.names(gspc_prices.returns)),
                y=-Apple_Garch2,
                color="Volatility"))+
  theme_bw()+
  theme1+
  ylab("AAPL log returns")+
  ggtitle("GARCH(1,2)")

garch_plot 


























































