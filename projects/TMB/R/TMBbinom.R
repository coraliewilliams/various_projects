### Simulate some data
set.seed(666)
# The true beta values
true.beta <- c(5, 1.5, -3)

# design matrix, with intercept and two covariates
X <- cbind(Int = rep(1, 100), X1 = rnorm(100), X2 = rnorm(100, mean = 2))

# create the linear mean
XB <- (X %*% true.beta) 

# get the inv-logit function

e =  rnorm(100)
prob <- exp(XB)/(1 + exp(XB))

Y = rbinom(100, 1, prob)

df = data.frame(Y, X)

glm(Y ~ X - 1, data = df, family = "binomial")

setwd("~/Code/TMB_Tutorials/")

library(TMB)

compile("CPPbinom.cpp")
dyn.load(dynlib("CPPbinom"))

obj <- MakeADFun(data = list(y = Y, X = X),
                 parameters = list(beta = rep(0, dim(X)[2])),
                 DLL = "CPPbinom",
                 silent = TRUE
)

res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")

sdreport(obj)
