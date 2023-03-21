//GLLVMs for Poisson distribution
#include<math.h>
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //declares all data and parameters used
  DATA_MATRIX(y);
  DATA_MATRIX(x);
  DATA_INTEGER(num_lv);
  PARAMETER_VECTOR(b0);
  PARAMETER_MATRIX(b);
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(loglam);
  PARAMETER_MATRIX(u); //latent variables, u, are treated as parameters
  
  vector<Type> lam_diag = exp(loglam); 
  int n = y.rows();
  int p = y.cols();
  //To create lambda as matrix upper triangle
  matrix<Type> newlam(num_lv,p);
  for (int j = 0; j < p; j++){
    for (int i = 0; i < num_lv; i++){
      if (j < i)
        newlam(i, j) = 0;
      else if(i == j)
        newlam(i, j) = lam_diag(j);
      else
        newlam(i, j) = lambda(i*p - (i + 1)*i/2 + (j - 1) - i   );
    }
  }
  
  matrix<Type> lam = u*newlam;
  
  //eta function b0 + x*b + u*lambda
  matrix<Type> eta(n,p);
  eta = x*b + lam;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      eta(i, j) = b0(j) + eta(i, j);
    }
  }
  
  Type nll = 0.0; // initial value of log-likelihood
  //latent variable is assumed to be from N(0,1)
  for (int j = 0; j < u.cols(); j++){
    for (int i = 0; i < n; i++) {
      nll -= dnorm(u(i,j), Type(0), Type(1), true);
    }
  }
  //likelihood poisson model with the log link function
  for (int j = 0; j < p ; j++){
    for (int i = 0 ; i < n; i++) {
      nll -= dpois(y(i,j), exp(eta(i,j)), true);
    }
  }
  
  REPORT(newlam);
  REPORT(lambda);
  REPORT(lam);
  REPORT(b0);
  REPORT(b);
  REPORT(u);
  
  return nll;
}