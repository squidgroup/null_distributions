data{
  int<lower=0> N;
  int<lower=0> N_ID;
  int<lower=0, upper=1> y[N];
  int<lower=0,upper=N_ID> ID[N];
}
parameters{
  real beta_0;
  vector[N_ID] ID_effects_scaled;
  real<lower=0> sigma_ID;
}
model {
  vector[N] mu = beta_0 + sigma_ID * ID_effects_scaled[ID];

  ID_effects_scaled ~ normal(0,1);

  beta_0 ~ normal(0,10);
  sigma_ID~ cauchy(0,2);

  y ~ bernoulli_logit(mu);
}
generated quantities {
  real sigma2_ID = sigma_ID^2;
}
