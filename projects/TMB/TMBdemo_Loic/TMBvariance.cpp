// LT 26/05/2016
// Toy example compute variance with TMB

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // x: the data
  DATA_VECTOR(x);
  
  // mean
  PARAMETER(m);

  // variance
  PARAMETER(sigma2);
  
  Type nLL=0;
  
  nLL -= dnorm(x, m, sqrt(sigma2), true).sum();
    
  return nLL;
}
