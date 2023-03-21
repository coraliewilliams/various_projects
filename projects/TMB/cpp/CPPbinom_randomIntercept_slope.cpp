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
  DATA_INTEGER(ngroups);        // number of levels in random effects
  
  // Parameters
  PARAMETER_VECTOR(Beta);         // Vector of beta values
  PARAMETER_ARRAY(u);             // Intercept for given random effect (/factor)
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(transformed_rho);     // parameter of correlation
  
  // Load namespace which contains the multivariate distributions
  using namespace density;
  /// define a matrix for the var-covar matrix for the multivariate normal
  matrix<Type> cov(k_size, k_size); 
  vector<Type> sd = exp(logsig1);   /// Take the exponential of the log std dev.
  Type rho = 2.0 / (1.0 + exp(-transformed_rho)) - 1.0;   /// To keep the correlation coef between -1, 1, use a shifted logistic form
  cov(0,0) = pow(sd[0],2);
  cov(1,1) = pow(sd[1],2);
  cov(0,1) = rho*sd[1]*sd[0];
  cov(1,0) = rho*sd[1]*sd[0];
  
  Type nll = 0.0;                 // initialize negative log likelihood
  
  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector
  vector<Type> ZU = X * u; // pre-calculate the design matrix times beta vector
  
  Type Size = 1;
  int k;                   // will act as a loop control variable between R and cpp
  
  for(int i = 0; i < Y.size(); i++){
    k = Factor(i) - 1;       // set the LCV to reflect the factor level of the observations
    nll -= dbinom_robust(Y(i), Size, log(exp(XB(i)  + u(k)  )), true);
  }
  
  // Component 2 - Random effects distribution
  MVNORM_t<Type> neg_log_density(cov);
  vector<Type> uj(k_size);
  
  for(int j = 0; j < ngroups; j++){
    uj = u.col(j);
    nll += neg_log_density(uj); // Process likelihood
  }
  
  return nll;
}