### Simulate some data

# The true beta values
true.beta <- c(5, 1.5, -3)

# design matrix, with intercept and two covariates
X <- cbind(Int = rep(1, 100), X1 = rnorm(100), X2 = rnorm(100, mean = 2))
# Add a 4 level factor to the data
dat <- as.data.frame(cbind(X, X3 = rep(1:4, length.out = dim(X)[1])))
dat$X3 <- factor(dat$X3)

# create the response with residuals
Y <- (X %*% true.beta) + rnorm(100)

# Add random amounts to the response within each factor level
dat$Y <- NULL

for (i in levels(dat$X3)) {
  dat$Y[dat$X3 == i] <- Y[dat$X3 == i] + rnorm(1)
}

library(lme4)
# # Fit the random intercept model using lmer from lme4
ri.mod = lmer(Y ~ X1 + X2 + (1|X3), data = dat,REML = FALSE)

setwd("~/Code")

library(TMB)
compile("lmer.cpp") 
dyn.load("lmer") 

# Create the objective function - setting 'u' to be random
ni=length(levels(dat$X3))
k=length(true.beta)
obj <- MakeADFun(data = list(X3 = as.numeric(dat$X3), Y = dat$Y, 
                             X = as.matrix(dat[ , c("Int", "X1", "X2")])),
                 parameters = list(Beta = rep(0, k), 
                                   u = rep(0, ni), 
                                   logsig1 = 0, logsig0 = 0),
                 random = "u",  ##u are random effects, will be integrated out w. Laplace
                 DLL = "lmer",
                 silent = TRUE
)

res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", 
             control = list(maxit = 1000), hessian = T)
sdreport(obj)
