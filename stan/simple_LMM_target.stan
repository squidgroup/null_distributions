data{
  int<lower=0> N;
  int<lower=0> N_ID;
  vector[N] y;
  int<lower=0,upper=N_ID> ID[N];
}
parameters{
  real beta_0;
  vector[N_ID] ID_effect_scaled;
  real<lower=0> sigma_ID;
  real<lower=0> sigma_E;
}
model{
  vector[N] mu = beta_0 + sigma_ID * ID_effects_scaled[ID];

  target += normal_lpdf( ID_scaled| 0, 1);
  target += normal_lpdf( beta_0| 0, 10);
  target += cauchy_lpdf( sigma_ID| 0, 2);
  target += cauchy_lpdf( sigma_E| 0, 2);
  target += normal_lpdf( y| mu, sigma_E);

}
