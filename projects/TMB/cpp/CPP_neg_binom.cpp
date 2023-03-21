// Simple Random Intercept Model
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data to be input
  DATA_VECTOR(Y);         // Response vector
  DATA_MATRIX(X);         // Design matrix
  DATA_MATRIX(Z);         // Random effect matrix
  DATA_IVECTOR(group);        // The factor for which we require random intercepts
  DATA_INTEGER(k_size);        // number of random effects
  DATA_INTEGER(nlevels);        // number of levels in random effects
  
  // Parameters
  PARAMETER_VECTOR(Beta);         // Vector of beta values
  PARAMETER_ARRAY(u);             // Intercept for given random effect (/factor)
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(logk);                // Dispersion parameter
  PARAMETER(transformed_rho);     // parameter of correlation
  
  // Load namespace which contains the multivariate distributions
  using namespace density;
  /// define a matrix for the var-covar matrix for the multivariate normal
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
  
  int N = Y.size();
  
  vector<Type> mu(N);
  vector<Type> eta(N);
  vector<Type> Zi(k_size);
  vector<Type> uj(k_size);
  int k;                   // will act as a loop control variable between R and cpp
  
  Type k_disp = exp(logk);
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector

  for(int i = 0; i < N; i++){
    k = group(i) - 1;       // set the LCV to reflect the group level of the observations
    uj = u.col(k);
    Zi = Z.row(i);
    // eta
    eta(i) = XB(i) + (Zi * uj).sum();
  }
  
  mu = exp(eta);
  
  // // Component 1 -  Observations: E(X|u)= nu(X|u)= XBeta + Zu
  Type nll = 0.0;                 // initialize negative log likelihood
  for(int i = 0; i < N; i++){
    nll -= dnbinom2(Y(i), mu(i), mu(i) + k_disp*pow(mu(i),2) , true);
  }
  
  // Component 2 - Random effects distribution
  MVNORM_t<Type> neg_log_density(covrand);
  for(int j = 0; j < nlevels; j++){
    uj = u.col(j);
    nll += neg_log_density(uj); // Process likelihood
  }
  
  ADREPORT(covrand);
  REPORT(covrand);
  ADREPORT(sd);
  REPORT(sd);
  ADREPORT(rho);
  REPORT(rho);
  ADREPORT(k_disp);
  REPORT(k_disp);
  
  return nll;
}