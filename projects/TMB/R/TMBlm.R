### Simulate some data

# The true beta values
true.beta <- c(5, 1.5, -3)

# design matrix, with intercept and two covariates
X <- cbind(Int = rep(1, 100), X1 = rnorm(100), X2 = rnorm(100, mean = 2))

# create the response with residuals
Y <- (X %*% true.beta) + rnorm(100)

lm(Y ~ X - 1)

setwd("~/Code/TMB_Tutorials")

library(TMB)

compile("CPPlm.cpp")
dyn.load(dynlib("CPPlm"))

obj <- MakeADFun(data = list(Y = Y, X = X),
                 parameters = list(Beta = rep(0, dim(X)[2]), logsig = 0),
                 DLL = "CPPlm",
                 silent = TRUE
)
obj$report



res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")

sdreport(obj)
