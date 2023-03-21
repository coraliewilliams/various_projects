// Simple Random Intercept Model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data to be input
  DATA_VECTOR(Y);         // Response vector
  DATA_MATRIX(X);         // Design matrix
  DATA_IVECTOR(Factor);        // The factor for which we require random intercepts
  DATA_INTEGER(k_size);        // number of random effects
  
  // Parameters
  PARAMETER_VECTOR(Beta);         // Vector of beta values
  PARAMETER_VECTOR(u);            // Intercept for given random effect (/factor)
  PARAMETER(logsig);      // Random effect sd
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(transformed_rho);     // parameter of correlation
  
  // Load namespace which contains the multivariate distributions
  using namespace density;

  matrix<Type> cov(k_size, k_size); /// define a matrix for the var-covar matrix for the multivariate normal
  vector<Type> sd = exp(logsig1);   /// Take the exponential of the log std dev.
  Type rho = 2.0 / (1.0 + exp(-transformed_rho)) - 1.0;   /// To keep the correlation coef between -1, 1, use a shifted logistic form

  for(int i = 0; i < k_size; i++){
    for(int j = 0; j < k_size; j++){
      if(i == j){
        cov(i, j) = sd[i];
      } else {
        cov(i, j) = rho*sd[i]*sd[j];
      }
    }
  }

  REPORT(cov)

  Type nll = 0.0;                 // initialize negative log likelihood
  // nll = MVNORM(cov)(u);
  
 // number of groups in the random effect
  // int ngroups = u.size();
  
  // // nll += MVNORM(cov)(u);
  // 
  // Type Size = 1; /// the number of bernoulli trials
  // int k;                   // will act as a loop control variable
  // 
  // // // Component 1 -  Observations: E(Y|u)= logit(Y|u)= XBeta + u
  // vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector
  // 
  // for(int i = 0; i < Y.size(); i++){
  //   k = Factor(i) - 1;       // set the LCV to reflect the factor level of the observations
  //   nll -= dbinom_robust(Y(i), Size, log(exp( XB(i)  ) ), true); //exp( XB(i)  + + u(k)
  // }

  int ngroups = u.size();  // define the number of factor levels i.e. random intercepts

  Type zero = 0.0;         // a constant
  int k;                   // will act as a loop control variable


  // Component 2 - Prior: intercept_j ~ N(0,sig1)
  for(int j = 0; j < ngroups; j++){
    nll -= dnorm(u(j), zero, exp(logsig), true);
  }

  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector

  Type Size = 1;

  for(int i = 0; i < Y.size(); i++){
    k = Factor(i) - 1;       // set the LCV to reflect the factor level of the observations
    nll -= dbinom_robust(Y(i), Size, log(exp(XB(i)  + u(k)  )), true);
  }

  return nll;
}
