//LT 25/05/2016
// A TMB version of negative binomial glmm
// for fast estimation of likelihood ratio null distribution

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // y: the response
  DATA_VECTOR(y);

  // X: design matrix of linear predictors
  DATA_MATRIX(X);

  // fixed effects parameters
  PARAMETER_VECTOR(beta);

  Type nLL = 0.0;
  
  vector<Type> XB = X * beta; // pre-calculate the design matrix times beta vector
  // vector<Type> mu = exp(XB)/(1 + exp(XB));
  
  Type Size = 1;

  for(int i=0; i<y.size(); i++){
    nLL -= dbinom_robust(y(i), Size, XB(i), true);
  }
    
  return nLL;
}
