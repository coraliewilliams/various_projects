# obj <- MakeADFun(data = data_mm,
#                  parameters = parameters,
#                  random = "b", ##u are random effects, will be integrated out w. Laplace
#                  DLL = "CPPGLLVM_poisson_long",
#                  silent = TRUE)
#
# opt <- optim(obj$par, obj$fn, obj$gr,
#              control = list(maxit = 10000), method = "BFGS")
#
# ### sdreport for estimates
# sd <- sdreport(obj)
#
# ### for maximum likelihood
# opt$value


# https://groups.google.com/forum/#!topic/tmb-users/ysLVa52rVsM
# If I return values from the TMB model using REPORT() as well as ADREPORT() the values differ by a very small percentage.
# the answer was that I wasn't explicitly calling the reported values using the best
# fitting parameters. My call to obj$report() should actually have been obj$report(obj$env$last.par.best).

################### ERRORS ###########
# Error in optim(obj$par, obj$fn, obj$gr, control = list(maxit = 10000),  :
#                  initial value in 'vmmin' is not finite
#  is due to poor initial starting parameter values and has a negative

#
# Error in optimHess(par.fixed, obj$fn, obj$gr) :
#   gradient in optim evaluated to length 1 not 204
