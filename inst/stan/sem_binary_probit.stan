data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  array[N, P] int<lower=0, upper=1> y;
  array[P, K] int<lower=0, upper=1> loading_allowed;
  matrix[P, K] loading_fixed;
  array[P, K] int<lower=0> loading_group;
  int<lower=0> G;
  real<lower=0> global_scale;
}
parameters {
  matrix[N, K] eta;
  matrix[P, K] loading_raw;
  vector[P] threshold;
  real<lower=0> tau;
  vector<lower=0>[G] kappa;
  matrix<lower=0>[P, K] psi;
}
transformed parameters {
  matrix[P, K] loading = loading_fixed;
  for (p in 1:P) for (k in 1:K) if (loading_allowed[p, k] == 1) {
    if (loading_group[p, k] > 0)
      loading[p, k] = loading_raw[p, k] * tau *
        kappa[loading_group[p, k]] * psi[p, k];
    else loading[p, k] = loading_raw[p, k];
  }
}
model {
  tau ~ normal(0, global_scale);
  kappa ~ cauchy(0, 1);
  to_vector(psi) ~ cauchy(0, 1);
  to_vector(loading_raw) ~ normal(0, 1);
  threshold ~ normal(0, 1.5);
  to_vector(eta) ~ normal(0, 1);
  for (n in 1:N) for (p in 1:P)
    y[n, p] ~ bernoulli(Phi(dot_product(loading[p], eta[n]) - threshold[p]));
}
generated quantities {
  array[N, P] int y_rep;
  for (n in 1:N) for (p in 1:P)
    y_rep[n, p] = bernoulli_rng(Phi(dot_product(loading[p], eta[n]) - threshold[p]));
}
