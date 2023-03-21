#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // The data to be input
  DATA_VECTOR(Y); // Response vector
  DATA_MATRIX(X); // Design matrix
  // The model parameters
  PARAMETER_VECTOR(Beta); // Vector of our 3 beta values
  PARAMETER(logsig); // natural log of the residual sd
  // require sigma > 0 therefore pass it to the objective transformation by the natural logarithm
  Type nll;
  nll = -sum(dnorm(Y, X*Beta, exp(logsig), true));
  // dnorm provides the gaussian pdf similar to R syntax
  // true is a logical that returns log of the density
  // We also take the negative of this value as our optimisation algorithms in R
  //   will default to minimisation (whereas we require maximisation).
  return nll;
  
}