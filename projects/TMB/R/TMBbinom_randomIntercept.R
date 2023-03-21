### Simulate some data
set.seed(666)

# The true beta values
true.beta <- c(5, 1.5, -3)

### number of observations
n.obs = 100

# design matrix, with intercept and two covariates
X <- cbind(Int = rep(1, n.obs), 
           X1 = rnorm(n.obs), 
           X2 = rnorm(n.obs, mean = 2))

### Adding in a covariate for a random factor
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
  dat$Y[dat$X3 == i] <- rbinom(size,1, prob)
}


library(lme4)
# # Fit the random intercept model using lmer from lme4
ri.mod = glmer(Y ~ X1 + X2 + (1|X3) , data = dat, family="binomial")


setwd("~/Code/TMB_Tutorials/")

library(TMB)

compile("CPPbinom_randomIntercept.cpp")
dyn.load(dynlib("CPPbinom_randomIntercept"))

# Create the objective function - setting 'u' to be random
k=length(true.beta)
ni=length(levels(dat$X3))

data = list(X3 = as.numeric(dat$X3), Y = dat$Y, X = as.matrix(dat[ , c("Int", "X1", "X2")]))
params = list(Beta = rep(0, k), u = rep(0, ni), logsig1 = 0)

obj <- MakeADFun(data = data,
                 parameters = params,
                 random = "u", ##u are random effects, will be integrated out w. Laplace
                 DLL = "CPPbinom_randomIntercept",
                 silent = TRUE)

res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", 
             control = list(maxit = 10000))

sdreport(obj)
ri.mod
