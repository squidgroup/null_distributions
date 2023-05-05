library(bridgesampling)
library(rstan)

### generate data ###
set.seed(12345)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n_pop=1

squid_dat <- simulate_population(
	data_structure= make_structure("ID(80)",repeat_obs=2),
	parameters= list(
		intercept=mu,
		ID=list(vcov=tau2), 
		residual=list(vcov=sigma2)),
	n_pop=n_pop,
)
data <- get_population_data(squid_dat)

# models
stancodeH0 <- 'data{
  int<lower=0> N;
  vector[N] y;
}
parameters{
  real beta_0;
  real<lower=0> sigma_E;
}
model{
  target += normal_lpdf( beta_0| 0, 10);
  target += cauchy_lpdf( sigma_E| 0, 2);
  target += normal_lpdf( y| beta_0, sigma_E);
}
'
stancodeH1 <- 'data{
  int<lower=0> N;
  int<lower=0> N_ID;
  vector[N] y;
  int<lower=0,upper=N_ID> ID[N];
}
parameters{
  real beta_0;
  vector[N_ID] ID_effects_scaled;
  real<lower=0> sigma_ID;
  real<lower=0> sigma_E;
}
model{
  vector[N] mu = beta_0 + sigma_ID * ID_effects_scaled[ID];

  target += normal_lpdf( ID_effects_scaled| 0, 1);
  target += normal_lpdf( beta_0| 0, 10);
  target += cauchy_lpdf( sigma_ID| 0, 2);
  target += cauchy_lpdf( sigma_E| 0, 2);
  target += normal_lpdf( y| mu, sigma_E);
}
'
# compile models
stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

stan_dat <- list(
	N = nrow(data),
	N_ID = length(unique(data[,"ID"])),
	y = data[,"y"],
	ID = data[,"ID"]
)

# fit models
stanfitH0 <- sampling(stanmodelH0, data = stan_dat,
                      iter = 15000, warmup = 2000, chains = 3, cores = 3)
stanfitH1 <- sampling(stanmodelH1, data = stan_dat,
                      iter = 15000, warmup = 2000, chains = 3, cores = 3)

# compute log marginal likelihood via bridge sampling for H0
H0.bridge <- bridge_sampler(stanfitH0, silent = TRUE)

# compute log marginal likelihood via bridge sampling for H1
H1.bridge <- bridge_sampler(stanfitH1, silent = TRUE)


# compute Bayes factor
BF01 <- bf(H1.bridge,H0.bridge)
print(BF01)

error_measures(H0.bridge)
error_measures(H1.bridge)

# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(H0.bridge, H1.bridge)
print(post1)