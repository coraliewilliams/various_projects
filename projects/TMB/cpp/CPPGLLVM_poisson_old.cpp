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
  PARAMETER_VECTOR(lambda2);
  PARAMETER_MATRIX(u); //latent variables, u, are treated as parameters
  
  int n = y.rows();
  int p = y.cols();
  //To create lambda as matrix upper triangle
  matrix<Type> newlam(num_lv,p);
  for (int j = 0; j < p; j++){
    for (int i = 0; i < num_lv; i++){
      if (j < i){
        newlam(i, j) = 0;
      } else{
        newlam(i, j) = lambda(j);
        if (i > 0){
          newlam(i, j) = lambda(i + j + i*p - (i*(i - 1)) / 2 - 2*i);
        }
      }
    }
  }
  
  matrix<Type> newlam2(num_lv,p);
  for (int j = 0; j < p; j++){
    for (int i = 0; i < num_lv; i++){
      if (j < i){
        newlam2(i, j) = 0;
      } else{
        newlam2(i, j) = lambda2(j);
        if (i > 0){
          newlam2(i, j) = lambda2(i + j + i*p - (i*(i - 1)) / 2 - 2*i);
        }
      }
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
  REPORT(newlam2);
  REPORT(lambda);
  REPORT(lam);
  REPORT(b0);
  REPORT(b);
  
  return nll;
}