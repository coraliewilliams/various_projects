setwd("~/Code/TMB_Tutorials/")

### Simulate some data
set.seed(666)

# The true beta values
true.beta <- c(5, 1.5, -3)

# number of observations
n.obs = 200

# design matrix, with intercept and two covariates
X <- cbind(Int = rep(1, n.obs), 
           X1 = rnorm(n.obs), 
           X2 = rnorm(n.obs, mean = 2))

# Adding in a covariate for a random factor
dat <- as.data.frame(cbind(X, X3 = rep(1:4, length.out = dim(X)[1])))

dat$X3 <- factor(dat$X3)

# create the linear mean from just the covariates (i.e., not including the random factor)
XB <- (X %*% true.beta) 

# Add random amounts to the response within each factor level
dat$Y <- NULL
for (i in levels(dat$X3)) {
  size = sum(dat$X3==i)
  ## random error for level specific
  u = rnorm(1)
  prob <- exp(XB[dat$X3 == i] + u)/(1 + exp(XB[dat$X3 == i] + u))
  ## using the probability generate random binomial observations
  dat$Y[dat$X3 == i] <- rbinom(size, 1, prob)
}

library(lme4)
# # Fit the random intercept model using lmer from lme4
ri.mod = glmer(Y ~ X1 + X2 + (1|X3) , data = dat, family="binomial")
### Model
### g(E(y)) = eta = XB + Zu
### Y is n * 1, where n is the number of obsverations
### p is the number of fixed effects
### u ~ N(0, sigma), sigma is the variance-covariance matrix

#### ------------------------------------------------------------------
#### TMB Part 
#### ------------------------------------------------------------------
library(TMB)
# library(TMBdebug)

compile("CPPbinom_randomIntercept_Covarmatrix.cpp")

# Create the objective function
### First define the parameters needed for the function
p = length(true.beta) # number of fixed effects
n.obs = length(dat$Y) # number of obs
k = ncol(as.matrix(dat$X3)) # number of random effects
li = length(levels(dat$X3)) # number of levels in random effect

## Data for the MakeADFun
data_mm = list(Y = dat$Y,
               X = as.matrix(dat[ , c("Int", "X1", "X2")]),
               Factor = as.numeric(dat$X3), 
               k_size = k)

## parameters to estimate
params = list(Beta = rep(0, p),
              u = rep(0, li),
              logsig1 = rep(0.5, k),
              transformed_rho = rep(0,1))

#### ------------------------------------------------------------------
#### MakeADFun 
#### ------------------------------------------------------------------
dyn.load(dynlib("CPPbinom_randomIntercept_Covarmatrix"))

obj <- MakeADFun(data = data_mm,
                 parameters = params,
                 random = "u", ##u are random effects, will be integrated out w. Laplace
                 DLL = "CPPbinom_randomIntercept_Covarmatrix",
                 silent = TRUE)


res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", 
             control = list(maxit = 10000)
             )

sdreport(obj)

obj$report()$cov

ri.mod
