#include /inst/stan/include/otherstanfunctions.stan
data {
  int<lower=0> T; //
  int<lower=0> K; //
  vector[T] y;
  matrix[T,T*K] big_s;
  vector[T*K] alpha;
  matrix[T*K,T*K] omega;
  matrix[T,T]  sigma;
  int distrib_method;
}

parameters {
  vector[T*K] mu;
}

transformed parameters{
  vector[T*K] beta = mod_rp_fun(mu, distrib_method);
  vector[T] Xbeta = big_s * beta ;
}

model {
  mu ~ multi_normal(alpha, omega);
  y ~ multi_normal(Xbeta, sigma);
}

generated quantities{
  real lupost = multi_normal_lpdf(mu|alpha, omega) + multi_normal_lpdf(y|Xbeta, sigma);
}

