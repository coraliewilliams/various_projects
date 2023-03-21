# LT 17/10/2015

# see whether TMB could speed up glmm fitting for manglmmm

library(TMB)

setwd("~/Code/TMB_Tutorials/TMBdemo_Loic/")

compile("TMBvariance.cpp")
dyn.load(dynlib("TMBvariance"))

x = rnorm(20)

f = MakeADFun(
  data=list(x=x),
  parameters=list(m=0, sigma2=1),
  #andom = "u",
  DLL="TMBvariance", silent=F)

# fit model
tmbfit = do.call(optim, f)

# get results
cat("tmb ", tmbfit$par,"\n")
cat("mean(x), mean((x-m)^2) ", mean(x), ", ", mean((x-mean(x))^2), "\n")

# profile likelihood
prof = tmbprofile(f, "sigma2")
confint(prof)
plot(prof)

# now let's integrate stuff.. REML
g = MakeADFun(
  data=list(x=x),
  parameters=list(m=0, sigma2=1),
  random = "m",
  DLL="TMBvariance", silent=F)

tmbfit.reml = do.call(optim, g)

# get results
cat("tmb ", tmbfit.reml$par,"\n")
cat("var(x) ", var(x), "\n")

# for some reason profile doesnt work ?
prof.reml = tmbprofile(g, "sigma2")

# do it manually
LLs =  sapply(prof[,1], g$fn)
# standardize likelihood
LLs = LLs - min(LLs) + min(prof[,2])
lines(prof[,1], LLs, col=2)
