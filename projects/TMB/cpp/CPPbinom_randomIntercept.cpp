// Simple Random Intercept Model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data to be input
  DATA_IVECTOR(X3);        // The factor for which we require random intercepts
  DATA_VECTOR(Y);          // Response vector
  DATA_MATRIX(X);          // Design matrix
  
  // Parameters
  PARAMETER_VECTOR(Beta);  // Vector of our 3 beta values
  PARAMETER_VECTOR(u);     // Intercept for given X3
  PARAMETER(logsig1);      // Random effect sd
  
  int ngroups = u.size();  // define the number of factor levels i.e. random intercepts
  Type nll = 0.0;         // initialize negative log likelihood
  
  Type zero = 0.0;         // a constant
  int k;                   // will act as a loop control variable
  
  // Component 2 - Prior: intercept_j ~ N(0,sig1)
  for(int j = 0; j < ngroups; j++){
    nll -= dnorm(u(j), zero, exp(logsig1), true);
  }

  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector

  Type Size = 1;

  for(int i = 0; i < Y.size(); i++){
    k = X3(i) - 1;       // set the LCV to reflect the factor level of the observations
    nll -= dbinom_robust(Y(i), Size, log(exp(XB(i)  + u(k)  )), true);
  }
  return nll;
}
