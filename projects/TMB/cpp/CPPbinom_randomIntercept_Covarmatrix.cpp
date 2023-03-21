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
  PARAMETER_VECTOR(u);             // Intercept for given random effect (/factor)
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(transformed_rho);     // parameter of correlation
  
  // Load namespace which contains the multivariate distributions
  using namespace density;
  /// define a matrix for the var-covar matrix for the multivariate normal
  matrix<Type> cov(k_size, k_size); 
  vector<Type> sd = exp(logsig1);   /// Take the exponential of the log std dev.
  Type rho = 2.0 / (1.0 + exp(-transformed_rho)) - 1.0;   /// To keep the correlation coef between -1, 1, use a shifted logistic form

  for(int i = 0; i < k_size; i++){
    for(int j = 0; j < k_size; j++){
      if(i == j){
        cov(i, j) = pow(sd[i], 2);
      } else {
        cov(i, j) = rho*sd[i]*sd[j];
      }
    }
  }
  
  Type nll = 0.0;                 // initialize negative log likelihood

  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector

  Type Size = 1;
  int k;                   // will act as a loop control variable between R and cpp
  
  // Working code
  for(int i = 0; i < Y.size(); i++){
    k = Factor(i) - 1;       // set the LCV to reflect the factor level of the observations
    nll -= dbinom_robust(Y(i), Size, log(exp(XB(i)  + u(k)  )), true);
  }

  // Component 2 - Random effects distribution
  MVNORM_t<Type> neg_log_density(cov);
  vector<Type> uj(1);

  int ngroups = u.size();  // define the number of factor levels i.e. random intercepts
  for(int j = 0; j < ngroups; j++){
    uj(0) = u(j);
    nll += neg_log_density(uj); // Process likelihood
  }

  
  return nll;
}
