data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  matrix[N, P] y;
  array[P, K] int<lower=0, upper=1> loading_allowed;
  matrix[P, K] loading_fixed;
  array[P, K] int<lower=0> loading_group;
  int<lower=0> G;
  array[K, K] int<lower=0, upper=1> structural_allowed;
  real<lower=0> global_scale;
}
parameters {
  matrix[N, K] eta;
  matrix[P, K] loading_raw;
  vector[P] intercept;
  vector<lower=0>[P] sigma;
  matrix[K, K] beta_raw;
  real<lower=0> tau;
  vector<lower=0>[G] kappa;
  matrix<lower=0>[P, K] psi;
}
transformed parameters {
  matrix[P, K] loading = loading_fixed;
  matrix[K, K] beta = rep_matrix(0, K, K);
  for (p in 1:P) for (k in 1:K) if (loading_allowed[p, k] == 1) {
    if (loading_group[p, k] > 0)
      loading[p, k] = loading_raw[p, k] * tau *
        kappa[loading_group[p, k]] * psi[p, k];
    else loading[p, k] = loading_raw[p, k];
  }
  if (K >= 2) for (k in 2:K) for (j in 1:(k - 1))
    if (structural_allowed[k, j] == 1) beta[k, j] = beta_raw[k, j];
}
model {
  tau ~ normal(0, global_scale);
  kappa ~ cauchy(0, 1);
  to_vector(psi) ~ cauchy(0, 1);
  to_vector(loading_raw) ~ normal(0, 1);
  to_vector(beta_raw) ~ normal(0, 0.5);
  intercept ~ normal(0, 2);
  sigma ~ normal(0, 1);
  for (n in 1:N) {
    eta[n, 1] ~ normal(0, 1);
    for (k in 2:K) eta[n, k] ~ normal(dot_product(beta[k, 1:(k - 1)],
                                                  eta[n, 1:(k - 1)]), 1);
    y[n] ~ normal(intercept + loading * to_vector(eta[n]), sigma);
  }
}
generated quantities {
  matrix[N, P] y_rep;
  for (n in 1:N) for (p in 1:P)
    y_rep[n, p] = normal_rng(intercept[p] + dot_product(loading[p], eta[n]),
                             sigma[p]);
}
