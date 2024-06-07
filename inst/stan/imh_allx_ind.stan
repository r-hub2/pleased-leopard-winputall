#include /inst/stan/include/otherstanfunctions.stan

data {
  int<lower=0> T; //
  int<lower=0> K; //
  int<lower=0> nb_sim; //
  vector[T] y;
  matrix[T,T*K] big_s;
  vector[T*K] alpha;
  matrix[T*K,T*K] omega;
  matrix[T,T]  sigma;
  int distrib_method;
  vector[T*K] start_value;
  vector[T*K] m_beta_importance;
  matrix[T*K, T*K] omega_importance;
}

model {

}

generated quantities{

  matrix[nb_sim, T*K] beta_accept;
  vector[nb_sim] log_p_mean;
  vector[T*K] beta_proposed = start_value;
  beta_accept[1] = beta_proposed' ;// start value

  real log_gg = multi_normal_lpdf(beta_proposed | m_beta_importance, omega_importance);

  vector[T*K] tbeta = mod_rp_fun(beta_proposed, distrib_method);
  vector[T] Xbeta = big_s * tbeta ;
  real log_pp_obs = multi_normal_lpdf(y | Xbeta, sigma);
  real log_pp_mix = multi_normal_lpdf(beta_proposed | alpha, omega) ;
  real log_pp = log_pp_obs + log_pp_mix;

  real p_t;
  real p_p;
  real R;
  real u;
  p_t = log_pp - log_gg;
  log_p_mean[1] = log_pp;

  int i = 2;
  while (i <= nb_sim) {
    beta_proposed = multi_normal_rng(m_beta_importance, omega_importance);
    log_gg = multi_normal_lpdf(beta_proposed | m_beta_importance, omega_importance);

    tbeta = mod_rp_fun(beta_proposed, distrib_method);
    Xbeta = big_s * tbeta ;

    log_pp_obs = multi_normal_lpdf(y | Xbeta, sigma);
    log_pp_mix = multi_normal_lpdf(beta_proposed | alpha, omega) ;

    log_pp = log_pp_obs + log_pp_mix;
    p_p = log_pp - log_gg;
    R   = fmin(1, exp(p_p - p_t));
    u   = uniform_rng(0, 1);
    if (u <= R) {
      // accept proposal
      beta_accept[i] = beta_proposed';
      p_t = p_p;  // save density
      log_p_mean[i]  = log_pp;
    } else {
      // stay with the current value
      beta_accept[i] = beta_accept[i - 1];
      log_p_mean[i]  = log_p_mean[i - 1];
      p_t = p_t * 1.0;
    }
    i += 1;
  }
}
