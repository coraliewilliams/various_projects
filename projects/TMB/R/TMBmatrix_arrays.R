### Simulate some data
set.seed(666)

true.beta <- c(5, 1.5, -3)

### number of observations
n.obs = 50

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
  ## individual error 
  # e = rnorm(size)
  prob <- exp(XB[dat$X3 == i] + u)/(1 + exp(XB[dat$X3 == i] + u))
  ## using the probability generate random binomial observations
  dat$Y[dat$X3 == i] <- rbinom(size, 1, prob)
}

dat$Yadd1 <- dat$Y + 1

# Create two vectors of different lengths.
vector1 <- c(5,9,3)
vector2 <- c(10,11,12,13,14,15)

# Take these vectors as input to the array.
arr1 <- array(c(vector1,vector2),dim = c(3,3,2))

setwd("~/Code/TMB_Tutorials/")

library(TMB)

compile("matrix_arrays.cpp")
dyn.load(dynlib("matrix_arrays"))

# Create the objective function - setting 'u' to be random
ni=length(levels(dat$X3))
k=length(true.beta)

obj <- MakeADFun(data = list(v1 = as.numeric(dat$Y)#, v2 = as.numeric(dat$Y), v3 = dat$Yadd1#,
                             
                             # m1 = as.matrix(dat[ , c("Int", "X1")]),
                             # m2 = as.matrix(dat[ , c("Int", "X2")]),
                             # m3 = as.matrix(dat[ , c("X1", "X2")]),
                             # m4 = as.matrix(dat[ , c("Int", "X1", "X2")]),
                             # a1 = arr1,
                             # a2 = arr1
                             ),
                 
                 parameters = list(list(0)),
                 
                 DLL = "matrix_arrays",
                 silent = TRUE
                 
)

res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", 
             control = list(maxit = 10000))

sdreport(obj)
ri.mod