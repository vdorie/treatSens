data {
  int<lower = 0> N;
  int<lower = 0> J;
  int<lower = 0> P;
  
  matrix[N, P] x;
  int<lower = 1, upper = J> g[N];
  int<lower = 0, upper = 1> z[N];
  vector[N] y;
  
  real zeta_z;
  real zeta_y;
  real<lower = 0.0, upper = 1.0> theta;
}
transformed data {
  vector<lower = 0.0, upper = 1.0>[N] z_real;
  
  z_real = to_vector(z);
}

parameters {
  vector[J] ranef_treatment;
  vector[J] ranef_response;
    
  real treatmentEffect;
  vector[P] beta_treatment;
  vector[P] beta_response;
  
  real<lower = 0.0> sigma_response;
  real<lower = 0.0> sigma_ranef_treatment;
  real<lower = 0.0> sigma_ranef_response;
}
transformed parameters {
  vector[N] log_p1;
  vector[N] log_p0;
  
  {
    real sigma_sq_response;
    vector[N] linpred_treatment;
    vector[N] propensity_1;
    vector[N] propensity_0;
    
    vector[N] mean_response;
    
    linpred_treatment = x * beta_treatment + ranef_treatment[g];
    
    for (n in 1:N) {
      propensity_1[n] = Phi(linpred_treatment[n] + zeta_z);
      propensity_0[n] = Phi(linpred_treatment[n]);
    }
    
    sigma_sq_response = sigma_response * sigma_response;
    mean_response = x * beta_response + ranef_response[g] + z_real * treatmentEffect;
    
    log_p1 = -0.5 * rows_dot_self(y - mean_response - zeta_y) / sigma_sq_response +
             log(z_real .* propensity_1 + (1.0 - z_real) .* (1.0 - propensity_1));
  
    log_p0 = -0.5 * rows_dot_self(y - mean_response) / sigma_sq_response +
             log(z_real .* propensity_0 + (1.0 - z_real) .* (1.0 - propensity_0));
  }
}
model {
  /* model is:
   * y | U, alpha ~ N(x * beta_y  + z * tau + u * zeta_y + G * alpha)
   * z | U, phi ~ Bern(pnorm(x * beta_z + u * zeta_z + G * phi)
   * U ~ Bern(pi)
   * alpha ~ N(0, sigma_alpha)
   * phi   ~ N(0, sigma_phi)
   *
   * However, since U is a discrete "parameter" which Stan cannot handle, we instead
   * integrate it out. This introduces a term into the joint posterior which looks like
   *   prod [ term{ u_i = 1 } + term{ u_i = 0 } ]
   * where the term contents are just the model as otherwise stated.
   *
   * For this we use the builtin function log_sum_exp(x, y) = log(exp(x) + exp(y)),
   * or more specifically, log_mix, which performs the same operations but uses our
   * mixing probability. Conseqently, log_p1 and log_p0 are the logs of the those terms
   * that we wish to sum.
   */
  
  for (n in 1:N)
    target += log_mix(theta, log_p1[n], log_p0[n]);
  target += -N * log(sigma_response);
  
  ranef_treatment ~ normal(0.0, sigma_ranef_treatment);
  ranef_response  ~ normal(0.0, sigma_ranef_response);
  
  treatmentEffect ~ student_t(3.0, 0.0, 4.0);
  beta_treatment ~ student_t(3.0, 0.0, 4.0);
  beta_response  ~ student_t(3.0, 0.0, 4.0);
  
  sigma_response ~ cauchy(0.0, 5.0);
  sigma_ranef_treatment ~ cauchy(0, 2.5);
  sigma_ranef_response ~ cauchy(0, 5.0);
}
generated quantities {
  vector[N] p;
  
  p = 1.0 ./ (1.0 + exp(log_p0 - log_p1));
}

