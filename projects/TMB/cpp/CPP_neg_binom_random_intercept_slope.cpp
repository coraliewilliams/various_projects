// Simple Random Intercept Model
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data to be input
  DATA_VECTOR(Y);         // Response vector
  DATA_MATRIX(X);         // Design matrix
  DATA_VECTOR(Z);         // Random effect matrix
  DATA_IVECTOR(Factor);        // The factor for which we require random intercepts
  DATA_INTEGER(k_size);        // number of random effects
  DATA_INTEGER(ngroups);        // number of levels in random effects
  
  // Parameters
  PARAMETER_VECTOR(Beta);         // Vector of beta values
  PARAMETER_ARRAY(u);             // Intercept for given random effect (/factor)
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(logk);         // Vector of dispersion parameters
  PARAMETER(transformed_rho);     // parameter of correlation
  
  // Load namespace which contains the multivariate distributions
  using namespace density;
  /// define a matrix for the var-covrandar matrix for the multivariate normal
  matrix<Type> covrand(k_size, k_size); 
  vector<Type> sd = exp(logsig1);   /// Take the exponential of the log std dev.
  Type rho = 2.0 / (1.0 + exp(-transformed_rho)) - 1.0;   /// To keep the correlation coef between -1, 1, use a shifted logistic form
  
  for(int i = 0; i < k_size; i++){
    for(int j = 0; j < k_size; j++){
      if(i == j){
        covrand(i, j) = pow(sd(i),2);
      } else {
        covrand(i, j) = rho*sd(i)*sd(j);
      }
    }
  }
  
  ADREPORT(covrand);
  REPORT(covrand);
  ADREPORT(sd);
  REPORT(sd);
  ADREPORT(rho);
  REPORT(rho);
  
  int N = Y.size();
  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector
  Type k_disp = exp(logk);
  vector<Type> b0(N);
  vector<Type> b1(N);
  vector<Type> mu(N);
  int k;                   // will act as a loop control variable between R and cpp
  
  for(int i = 0; i < N; i++){
    k = Factor(i) - 1;       // set the LCV to reflect the factor level of the observations
    b0(i) = u(0, k);
    b1(i) = u(1, k);
    // mu = exp(XB+ Zu)
    mu(i) = exp(XB(i)  + b0(i) + b1(i)*Z(i));
  }
  
  Type nll = 0.0;                 // initialize negative log likelihood
  for(int i = 0; i < N; i++){
    nll -= dnbinom2(Y(i), mu(i), mu(i) + k_disp*pow(mu(i),2) , true);
  }
  
  // Component 2 - Random effects distribution
  MVNORM_t<Type> neg_log_density(covrand);
  
  vector<Type> uj(k_size);
  for(int j = 0; j < ngroups; j++){
    uj = u.col(j);
    nll += neg_log_density(uj); // Process likelihood
  }
  return nll;
  
  }
//-------------------------------------------------------------