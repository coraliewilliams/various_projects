// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_data );
  DATA_INTEGER( n_factors );
  DATA_FACTOR( Factor );
  DATA_VECTOR( Y );     
  DATA_INTEGER( k_size );        // number of random effects
  
  // Parameters
  PARAMETER( X0 );
  PARAMETER( log_SD0 );
  PARAMETER_VECTOR(log_SDZ);      // Random effect sd
  PARAMETER_VECTOR( Z );          // Random effect
  // PARAMETER( transformed_rho );     // parameter of correlation
  
  // matrix<Type> cov(k_size, k_size); /// define a matrix for the var-covar matrix for the multivariate normal
  // vector<Type> sd = exp(log_SDZ);   /// Take the exponential of the log std dev.
  // Type rho = 2.0 / (1.0 + exp(-transformed_rho)) - 1.0;   /// To keep the correlation coef between -1, 1, use a shifted logistic form
  // 
  // for(int i = 0; i < k_size; i++){
  //   for(int j = 0; j < k_size; j++){
  //     if(i == j){
  //       cov(i, j) = sd[i];
  //     } else {
  //       cov(i, j) = rho*sd[i]*sd[j];
  //     }
  //   }
  // }
  
  // Objective funcction
  Type jnll = 0;
  
  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm( Y(i), X0 + Z(Factor(i)), exp(log_SD0), true );
  }

  // Probability of random coefficients
  for( int i=0; i<n_factors; i++){
    jnll -= dnorm( Z(i), Type(0.0), exp(log_SDZ[0]), true );
  }
  
  // using namespace density;
  // MVNORM_t<Type> neg_log_density(cov);
  // 
  //   jnll += neg_log_density(Z); // random effects likelihood

  
  // Reporting
  // Type SDZ = exp(log_SDZ);
  // Type SD0 = exp(log_SD0);
  // ADREPORT( SDZ );
  // REPORT( SDZ );
  // ADREPORT( SD0 );
  // REPORT( SD0 );
  // ADREPORT( Z );
  // REPORT( Z );
  // ADREPORT( X0 );
  // REPORT( X0 );
  // 
  // // bias-correction testing
  // Type MeanZ = Z.sum() / Z.size();
  // Type SampleVarZ = ( (Z-MeanZ) * (Z-MeanZ) ).sum();
  // Type SampleSDZ = pow( SampleVarZ + 1e-20, 0.5);
  // REPORT( SampleVarZ );
  // REPORT( SampleSDZ );  
  // ADREPORT( SampleVarZ );
  // ADREPORT( SampleSDZ );  
  
  return jnll;
}