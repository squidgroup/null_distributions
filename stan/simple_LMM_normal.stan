data{
  int<lower=0> N;
  int<lower=0> N_ID;
  vector[N] y;
  int<lower=0,upper=N_ID> ID[N];
  real<lower=0> normal_scale;
}
parameters{
  real beta_0;
  vector[N_ID] ID_effects_scaled;
  real<lower=0> sigma_ID;
  real<lower=0> sigma_E;
}
model{
  vector[N] mu = beta_0 + sigma_ID * ID_effects_scaled[ID];

  ID_effects_scaled ~ normal(0,1);

  beta_0 ~ normal(0,10);
  sigma_ID ~ normal(0,normal_scale);
  sigma_E ~ normal(0,normal_scale);

  y ~ normal(mu, sigma_E);
}
generated quantities{
  real sigma2_ID = sigma_ID^2;
  real sigma2_E = sigma_E^2;
}
