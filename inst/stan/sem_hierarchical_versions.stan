data {
  int<lower=1> N;
  int<lower=2> R;
  int<lower=2> K;
  int<lower=2> P;
  array[N] int<lower=1, upper=R> version;
  matrix[N, P] y;
  array[P, K] int<lower=0, upper=1> loading_allowed;
  matrix[P, K] loading_fixed;
  array[K, K] int<lower=0, upper=1> structural_allowed;
  int<lower=1, upper=3> pooling_mode; // 1=no, 2=complete, 3=partial
}
parameters {
  matrix[N, K] eta;
  matrix[P, K] loading_raw;
  vector[P] intercept;
  vector<lower=0>[P] sigma_y;
  matrix[K, K] beta_mu;
  matrix<lower=0>[K, K] beta_tau;
  array[R] matrix[K, K] beta_z;
}
transformed parameters {
  matrix[P, K] loading = loading_fixed;
  array[R] matrix[K, K] beta;
  for (p in 1:P) for (k in 1:K)
    if (loading_allowed[p, k] == 1) loading[p, k] = loading_raw[p, k];
  for (r in 1:R) {
    beta[r] = rep_matrix(0, K, K);
    for (k in 2:K) for (j in 1:(k - 1)) if (structural_allowed[k, j] == 1) {
      if (pooling_mode == 1) beta[r, k, j] = 0.5 * beta_z[r, k, j];
      else if (pooling_mode == 2) beta[r, k, j] = beta_mu[k, j];
      else beta[r, k, j] = beta_mu[k, j] + beta_tau[k, j] * beta_z[r, k, j];
    }
  }
}
model {
  to_vector(loading_raw) ~ normal(0, 1);
  intercept ~ normal(0, 2);
  sigma_y ~ normal(0, 1);
  to_vector(beta_mu) ~ normal(0, 0.5);
  to_vector(beta_tau) ~ normal(0, 0.25);
  for (r in 1:R) to_vector(beta_z[r]) ~ normal(0, 1);
  for (n in 1:N) {
    eta[n, 1] ~ normal(0, 1);
    for (k in 2:K)
      eta[n, k] ~ normal(dot_product(beta[version[n], k, 1:(k - 1)],
                                     eta[n, 1:(k - 1)]), 1);
    y[n] ~ normal(intercept + loading * to_vector(eta[n]), sigma_y);
  }
}
generated quantities {
  matrix[N, P] y_rep;
  for (n in 1:N) for (p in 1:P)
    y_rep[n, p] = normal_rng(intercept[p] + dot_product(loading[p], eta[n]),
                             sigma_y[p]);
}
