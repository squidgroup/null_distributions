data{
  int<lower=0> N;
  int<lower=0> N_ID;
  vector[N] y;
  vector[N] x_1;
  int<lower=0,upper=N_ID> ID[N];
}
parameters{
  real beta_0;
  real beta_1;
  vector[N_ID] int_effects_scaled;
  vector[N_ID] slope_effects_scaled;
  real<lower=0> sigma_int;
  real<lower=0> sigma_slope;
  real<lower=0> sigma_E;
}
model{
  vector[N] mu = beta_0 + beta_1*x_1 + sigma_int * int_effects_scaled[ID] + sigma_slope * slope_effects_scaled[ID] .*x_1;

  int_effects_scaled ~ normal(0,1);
  slope_effects_scaled ~ normal(0,1);

  beta_0 ~ normal(0,10);
  beta_1 ~ normal(0,10);
  sigma_int ~ cauchy(0,2);
  sigma_slope ~ cauchy(0,2);
  sigma_E ~ cauchy(0,2);

  y ~ normal(mu, sigma_E);
}
generated quantities{
  real sigma2_int = sigma_int^2;
  real sigma2_slope = sigma_slope^2;
}
