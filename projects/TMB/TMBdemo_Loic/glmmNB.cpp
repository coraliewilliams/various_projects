//LT 25/05/2016
// A TMB version of negative binomial glmm
// for fast estimation of likelihood ratio null distribution

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // y: the response
  DATA_VECTOR(y);

  // X: design matrix of linear predictors
  DATA_MATRIX(X);

  // Z: design matrix of ranefs
  DATA_SPARSE_MATRIX(Z);

  // sparse cholesky factor
  DATA_SPARSE_MATRIX(Lambda);

  // indicators for variance components
  DATA_IVECTOR(Lind);

  // variance components parameters
  PARAMETER_VECTOR(theta);

  // fixed effects parameters
  PARAMETER_VECTOR(beta);

  // conditional mode of the random effects
  PARAMETER_VECTOR(u);

  // overdispersion parameter of the negative binomial
  PARAMETER(alpha);

  Type nLL=0;

  //set the variance components at the right place in lambda
  // assigning point to indicators for variance components
  int    *lipt = Lind.data();
  // 
  Type *LamX = Lambda.valuePtr(), *thpt = theta.data();
  for (int i = 0; i < Lind.size(); ++i) {
    LamX[i] = thpt[lipt[i] - 1];
  }

  // contribution of ranefs to likelihood
  nLL -= dnorm(u, Type(0), Type(1), true).sum();

  // mu = exp(eta)
  vector<Type> mu = exp(X*beta+Z*(Lambda*u));

  // some debug
  //std::cout << "mu: " << mu << std::endl;

  vector<Type> variances = mu + mu.cwiseProduct(mu)/alpha;

  // parametrized with mean and variance
  nLL -= dnbinom2(y, mu, variances, true).sum();

  return nLL;
}
