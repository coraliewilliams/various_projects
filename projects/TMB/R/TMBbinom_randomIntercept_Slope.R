### Simulate some data
set.seed(666)

#First generate data at the unit level, ai ~ N(0, 1)
N <- 30
unit.df = data.frame(unit = c(1:N), a = rnorm(N))

head(unit.df, 3)

#Every unit will have its own slope, alpha_i, and intercept, Beta_i. The key is, these are now linearly related to ai. 

unit.df <-  within(unit.df, {
  E.alpha.given.a <-  1 - 0.15*a
  E.beta.given.a <-  3 + 0.3*a} )
head(unit.df, 3)

library(mvtnorm)
q = 0.2
r = 0.2
s = 0.5
cov.matrix <-  matrix(c(q^2, 
                        r*q*s, r*q*s,
                        s^2), nrow = 2, byrow = TRUE)
random.effects <-  rmvnorm(N, mean = c(0, 0), sigma = cov.matrix)

# Based on our choices of q and s above, slopes vary more than intercepts. 
# Also we've specified a correlation of .9 (i.e., .09/sqrt(.04???.25))
# so when the intercept's random effect is above its mean (of zero), the same is likely true for the slope. 
unit.df$alpha = unit.df$E.alpha.given.a + random.effects[, 1]
unit.df$beta = unit.df$E.beta.given.a + random.effects[, 2]

# We now move to the within-unit level, generating a fixed x-grid to be shared among units. 
J = 30
M = J * N  #Total number of observations
x.grid = seq(-4, 4, by = 8/J)[0:30]

within.unit.df <-  data.frame(unit = sort(rep(c(1:N), J)), j = rep(c(1:J), N), x =rep(x.grid, N))
head(within.unit.df)
# Next, we generate data from the individual linear models based on Alpha_i, Beta_i, and within-unit level noise ??ij???N(0,.752). 
flat.df = merge(unit.df, within.unit.df)

flat.df  <- within(flat.df, prob <- exp(x*beta)/(1 + exp(x*beta)))

## using the probability generate random binomial observations
flat.df  <- within(flat.df, y <- rbinom(M, 1, prob))

simple.df <-  flat.df[, c("unit", "a", "x", "y")]
head(simple.df, 3)
simple.df <- cbind(Int = rep(1, M), simple.df)

library(lme4)
my.lmer <-  glmer(y ~ x + (1 + x | unit), data = simple.df, family = "binomial")

dat <- simple.df
# > names(dat)
# [1] "Int"  "unit" "a"    "x"    "y"  
setwd("~/Code/TMB_Tutorials/")

library(TMB)

compile("CPPbinom_randomIntercept_slope.cpp")
dyn.load(dynlib("CPPbinom_randomIntercept_slope"))

# Create the objective function
### First define the parameters needed for the function
p = ncol(as.matrix(dat[ , c("Int", "x")])) # number of fixed effects
n.obs = length(dat$y) # number of obs
k = ncol(as.matrix(cbind(dat$Int, dat$x))) # number of random effects
li = length(unique(dat$unit)) # number of levels in random effect

## Data for the MakeADFun
data_mm = list(Y = dat$y,
            X = as.matrix(dat[ , c("Int", "x")]),
            Factor = as.numeric(dat$unit), 
            k_size = k,
            # k_size = 1, ## for when there is only random intercept
            ngroups = li)

## parameters to estimate
params = list(Beta = rep(0, p),
              u = array(rep(0, li), dim = c(2, li)),
              logsig1 = rep(0, k),
              # u = rep(0, li), ## for when there is only random intercept
              # logsig1 = rep(0, 1), ## for when there is only random intercept
              transformed_rho = rep(0,1)
              )

#### ------------------------------------------------------------------
#### MakeADFun 
#### ------------------------------------------------------------------
obj <- MakeADFun(data = data_mm,
                 parameters = params,
                 random = "u", ##u are random effects, will be integrated out w. Laplace
                 DLL = "CPPbinom_randomIntercept_slope",
                 silent = TRUE)

res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, #method = "BFGS", 
             control = list(maxit = 10000))

sdreport(obj)

