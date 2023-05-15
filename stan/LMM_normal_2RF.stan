data{
  int<lower=0> N;
  int<lower=0> N_group1;
  int<lower=0> N_group2;
  vector[N] y;
  int<lower=0,upper=N_group1> group1[N];
  int<lower=0,upper=N_group2> group2[N];

  real<lower=0> cauchy_scale;
}
parameters{
  real beta_0;
  vector[N_group1] group1_effects_scaled;
  vector[N_group2] group2_effects_scaled;
  real<lower=0> sigma_group1;
  real<lower=0> sigma_group2;
  real<lower=0> sigma_E;
}
model{
  vector[N] mu = beta_0 + sigma_group1 * group1_effects_scaled[group1] + sigma_group2 * group2_effects_scaled[group2];

  group1_effects_scaled ~ normal(0,1);
  group2_effects_scaled ~ normal(0,1);

  beta_0 ~ normal(0,10);
  sigma_group1 ~ cauchy(0,cauchy_scale);
  sigma_group2 ~ cauchy(0,cauchy_scale);
  sigma_E ~ cauchy(0,cauchy_scale);

  y ~ normal(mu, sigma_E);
}
generated quantities{
  real sigma2_group1 = sigma_group1^2;
  real sigma2_group2 = sigma_group2^2;
  real sigma2_E = sigma_E^2;
}
