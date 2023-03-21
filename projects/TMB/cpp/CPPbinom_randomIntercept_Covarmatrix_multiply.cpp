/ Simple Random Intercept Model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data to be input
  DATA_IVECTOR(X3);        // The factor for which we require random intercepts
  DATA_VECTOR(Y);          // Response vector
  DATA_MATRIX(X);          // Design matrix
  DATA_MATRIX(Z);         // Matrix for the random effect vars
  
  // Parameters
  PARAMETER_VECTOR(Beta);  // Vector of our 3 beta values
  PARAMETER_VECTOR(u);     // Intercept for given X3
  PARAMETER_VECTOR(logsig1);      // Random effect sd
  PARAMETER(transformed_rho); // parameter of correlation
  
  int k_size = Z.cols();       // define the number of random effects
  Type nll = 0.0;         // initialize negative log likelihood
  
  using namespace density;
  // 
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
  
  nll = MVNORM(cov)(u);
  
  Type Size = 1; /// the number of bernoulli trials
  int k;                   // will act as a loop control variable

  // // Component 1 -  Observations: E(X|u)= logit(X|u)= XBeta + u
  vector<Type> XB = X * Beta; // pre-calculate the design matrix times beta vector

  for(int i = 0; i < Y.size(); i++){
    k = X3(i) - 1;       // set the LCV to reflect the factor level of the observations
    nll -= dbinom_robust(Y(i), Size, log(exp(XB(i)  + Z(i)*u(k)  )), true);
  }

  return nll;
}
