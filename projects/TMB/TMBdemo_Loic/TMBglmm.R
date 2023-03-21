# LT 24/05/2016

# test whether TMB significantly speeds up newglmm simulation of LR null distribution

library(TMB)
# setwd("H:/David_Emails/Code/TMBdemo_Loic/")
compile("glmmNB.cpp")
dyn.load(dynlib("glmmNB"))

library(mvabund)
data(Tasmania)
library(lme4)

y= as.vector(Tasmania$copepods)
X = data.frame(spp=as.vector(col(Tasmania$copepods)), block=rep(Tasmania$block, 12),
               treat = rep(Tasmania$treatment, 12))

# fit null and alternative model (test for treatment effect)
m.0 = glmer.nb(y ~ 1 + (1|spp) + (1|block:spp), data=X)
m.1 = glmer.nb(y ~ 1 + treat + (treat|spp) + (1|block:spp), data=X)

# null model
f = MakeADFun(
  data=  getME(m.0, c("y","X","Z","Lambda", "Lind")),
  parameters=c(getME(m.0,c("theta","beta", "u")), alpha=1),
  random = "u",
  DLL="glmmNB", silent=F)

tmbfit.0 = nlminb(f$par,f$fn,f$gr, control=list(eval.max=10000, iter.max=5000))

# alternative
g = MakeADFun(
  data=  getME(m.1, c("y","X","Z","Lambda", "Lind")),
  parameters=c(getME(m.1,c("theta","beta", "u")), alpha=1),
  random = "u",
  DLL="glmmNB", silent=F)

tmbfit.1 = nlminb(g$par,g$fn,g$gr, control=list(eval.max=10000, iter.max=5000))

# profile on dispersion parameter ?
tmbprof = tmbprofile(g, "alpha")
plot(tmbprof)

# we could for example fit m.0 and m.h1 with the same alpha (impossible with glmer.nb)
# and test versus model with different alpha

# speed test a null distribution

simdata = simulate(m.0,10, use.u=T)

system.time({
  nd.tmb = apply(simdata, 2, function(ysim) {

    f = MakeADFun(
      data=  c(list(y=ysim), getME(m.0, c("X","Z","Lambda", "Lind"))),
      parameters=c(getME(m.0,c("theta","beta", "u")), alpha=1),
      random = "u",
      DLL="glmmNB", silent=T)
  
  tmbfit.0 = nlminb(f$par,f$fn,f$gr, control=list(eval.max=10000, iter.max=5000))
  
  # alternative
  g = MakeADFun(
    data=  c(list(y=ysim), getME(m.1, c("X","Z","Lambda", "Lind"))),
    parameters=c(getME(m.1,c("theta","beta", "u")), alpha=1),
    random = "u",
    DLL="glmmNB", silent=T)
  
  tmbfit.1 = nlminb(g$par,g$fn,g$gr, control=list(eval.max=10000, iter.max=5000))

  tmbfit.0$obj - tmbfit.1$obj
})
})
# 70s for 100 simulations

system.time({
nd.glmer = apply(simdata, 2, function(ysim) {
  m.0 = glmer.nb(ysim ~ 1 + (1|spp) + (1|block:spp), data=X)
  m.1 = glmer.nb(ysim ~ 1 + treat + (treat|spp) + (1|block:spp), data=X)
  logLik(m.1) - logLik(m.0)
})
})
# 731 s,  more than 50 warnings

# compare distributions
plot(nd.glmer, nd.tmb)
abline(a=0,b=1,col=2)

# try refit instead
system.time({
nd.glmer.fast = apply(simdata, 2, function(ysim) {
  m.0 = refit(m.0, ysim)
  m.1 = refit(m.1, ysim)
  logLik(m.1) - logLik(m.0)
})
})
# 80 s, 13 warnings

hist(nd.glmer.fast)

#only way to get something decent is to refit with same dispersion parameter
m.1=update(m.1, family=negative.binomial(getME(m.0, "glmer.nb.theta")))
system.time({
  nd.glmer.fast2 = apply(simdata, 2, function(ysim) {
    m.0 = refit(m.0, ysim)
    m.1 = refit(m.1, ysim)
    logLik(m.1) - logLik(m.0)
  })
})

# 81 s, 10 warnings
plot(nd.glmer.fast2, nd.glmer)
abline(a=0,b=1, col=2)
# yuck ! makes the test conservative




