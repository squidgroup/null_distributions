rm(list=ls())

library(parallel)
library(squidSim)
library(arm)
library(rstan)
rstan_options(auto_write = TRUE)

# setwd("/home/yi/Desktop/Perm_RRExample")
setwd("/Users/joelpick/github/null_distributions")
source("R/RR_functions.R")

#Simulated example
##Generate "real" simulated data set
set.seed("20230522")

squid_data <- simulate_population(
  data_structure=make_structure("individual(300)",repeat_obs=4),
  parameters = list(
    individual = list(
      names = c("ind_int","ind_slope"), 
      beta = c(1,0),
      vcov = c(0.2,0.1)
    ),
    observation= list(
      names = c("environment"),
      beta = c(0.2)
    ), 
    residual = list(
      vcov = c(0.7)
    ),
    interactions = list(
      names = c("ind_slope:environment"),
      beta = c(1)
    )
  )
)


data <- get_population_data(squid_data)

##Analyze real data set
stan_data = list(n = nrow(data),
                 n_ind=max(data$individual),
                 individual=data$individual,
                 x=data$environment,
                 y=data$y)

params<-c("beta_0", "beta_1", "sigma2_I1", "sigma2_I2", "sigma2_e")
stan_md <- stan_model(file = "Stan/RR_example.stan")

obs_md_sim <- sampling(stan_md, data = stan_data, pars = params, chains = 3, iter = 5500,  warmup = 500, thin = 5, cores=3)
obs_est_sim<-sapply(extract(obs_md_sim)[params],median)

##Permutation tests
n_perm<-100
t1<-system.time({  
  l_rand_ID_data<-ran_ID(data, n_perm)
  res_ID<-unlist(mclapply(1:n_perm,function(x){
    cat(x," ")
    mod(l_rand_ID_data[[x]])
  } , mc.cores=8))
})

l_rand_y_data<-ran_y(data, n_perm)
res_y<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(l_rand_y_data[[x]])
} , mc.cores=8))

l_rand_x_data<-ran_x(data, n_perm)
res_x<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(l_rand_x_data[[x]])
} , mc.cores=8))

l_rand_x_ID_data<-ran_x_ID(data, n_perm)
res_x_ID<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(l_rand_x_ID_data[[x]])
} , mc.cores=8))

##Bootstrapping - with no slope variance
t2<-system.time({  
squid_data_null_ID <- simulate_population(
  data_structure=make_structure("individual(300)",repeat_obs=4), 
  n_pop = n_perm,
  parameters = list(
    intercept = obs_est_sim[1],
    individual = list(
      names = c("ind_int","ind_slope"), 
      beta = c(1,0),
      vcov = c(obs_est_sim[3],0)
    ),
    residual = list(
      vcov = c(obs_est_sim[4]+obs_est_sim[5])
    ),
    interactions = list(
      names = c("ind_slope:environment"),
      beta = c(1)
    )
  ),
  known_predictors = list(
    predictors = data.frame(environment=data[,"environment"]), 
    beta = obs_est_sim[2]
  )
)

data_null_ID <- get_population_data(squid_data_null_ID, list=TRUE)

res_boot_ID<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(data_null_ID[[x]])
} , mc.cores=8))

})




##Bootstrapping - with no intercept or slope variance
squid_data_null_both <- simulate_population(
  data_structure=make_structure("individual(300)",repeat_obs=4), 
  n_pop = n_perm,
  parameters = list(
    intercept = obs_est_sim[1],
    individual = list(
     names = c("ind_int","ind_slope"), 
     beta = c(1,0),
     vcov = c(0,0)
    ),
    residual = list(
      vcov = c(obs_est_sim[3]+obs_est_sim[4]+obs_est_sim[5])
    ),
    interactions = list(
      names = c("ind_slope:environment"),
      beta = c(1)
    )
  ),
  known_predictors = list(
    predictors = data.frame(environment=data[,"environment"]), 
    beta = obs_est_sim[2]
  )
)

data_null_both <- get_population_data(squid_data_null_both, list=TRUE)

res_boot_both<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(data_null_both[[x]])
} , mc.cores=8))



##Results data frame
res<-data.frame(
  est=c(res_y, res_ID, res_x, res_x_ID, res_boot_ID, res_boot_both), 
  type=as.factor(rep(c("perm_y", "perm_ID", "perm_x", "perm_x_ID","boot_ID","boot_both"), each=n_perm)))


##P-values
tapply(res$est,res$type,function(x)mean(x>obs_est_sim["sigma2_I2"]))



### ---------------
##Real example
### ---------------
set.seed("20230522")

data_GT<-read.csv("Data/Raw/data_GT.csv")
stan_data_GT = list(n = nrow(data_GT),
                 n_ind=max(data_GT$individual),
                 individual=data_GT$individual,
                 x=data_GT$environment,
                 y=data_GT$y)

params<-c("beta_0", "beta_1", "sigma2_I1", "sigma2_I2", "sigma2_e")
obs_md_GT <- sampling(stan_md, data = stan_data_GT, pars = params, chains = 3, iter = 5500,  warmup = 500, thin = 30, cores=3)
obs_est_real<-sapply(extract(obs_md_GT)[params],median)

##Permutation tests
t3<-system.time({  
l_rand_ID_data_GT<-ran_ID(data_GT, n_perm)
res_ID_GT<-unlist(mclapply(1:n_perm,function(x){
    cat(x," ")
    mod(l_rand_ID_data_GT[[x]])
  } , mc.cores=8))
})

l_rand_y_data_GT<-ran_y(data_GT, n_perm)
res_y_GT<-unlist(mclapply(1:n_perm,function(x){
    cat(x," ")
    mod(l_rand_y_data_GT[[x]])
  } , mc.cores=8))

l_rand_x_data_GT<-ran_x(data_GT, n_perm)
res_x_GT<-unlist(mclapply(1:n_perm,function(x){
    cat(x," ")
    mod(l_rand_x_data_GT[[x]])
  } , mc.cores=8))

l_rand_x_ID_data_GT<-ran_x_ID(data_GT, n_perm)
res_x_ID_GT<-unlist(mclapply(1:n_perm,function(x){
    cat(x," ")
    mod(l_rand_x_ID_data_GT[[x]])
  } , mc.cores=8))

##Bootsraping
t4<-system.time({  
squid_data_null_ID_GT <- simulate_population(
  data_structure=data.frame(individual=data_GT$individual), 
  n_pop = n_perm,
  parameters = list(
    intercept = obs_est_real[1],
    individual = list(
      names = c("ind_int","ind_slope"), 
      beta = c(1,0),
      vcov = c(obs_est_real[3],0)
    ),
    residual = list(
      vcov = c(obs_est_real[4]+obs_est_real[5])
    ),
    interactions = list(
      names = c("ind_slope:environment"),
      beta = c(1)
    )
  ),
  known_predictors = list(
    predictors = data.frame(environment=data_GT[,"environment"]), 
    beta = obs_est_real[2]
  )
)

data_null_ID_GT <- get_population_data(squid_data_null_ID_GT, list=TRUE)

res_boot_ID_GT<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(data_null_ID_GT[[x]])
} , mc.cores=8))
})

squid_data_null_both_GT <- simulate_population(
  data_structure=data.frame(individual=data_GT$individual), 
  n_pop = n_perm,
  parameters = list(
    intercept = obs_est_real[1],
    individual = list(
      names = c("ind_int","ind_slope"), 
      beta = c(1,0),
      vcov = c(0,0)
    ),
    residual = list(
      vcov = c(obs_est_real[3]+obs_est_real[4]+obs_est_real[5])
    ),
    interactions = list(
      names = c("ind_slope:environment"),
      beta = c(1)
    )
  ),
  known_predictors = list(
    predictors = data.frame(environment=data_GT[,"environment"]), 
    beta = obs_est_real[2]
  )
)

data_null_both_GT <- get_population_data(squid_data_null_both_GT, list=TRUE)

res_boot_both_GT<-unlist(mclapply(1:n_perm,function(x){
  cat(x," ")
  mod(data_null_both_GT[[x]])
} , mc.cores=8))

res_GT<-data.frame(
  est=c(res_y_GT, res_ID_GT, res_x_GT, res_x_ID_GT, res_boot_ID_GT, res_boot_both_GT), 
  type=as.factor(rep(c("perm_y", "perm_ID", "perm_x", "perm_x_ID","boot_ID","boot_both"), each=n_perm)))


##P-values
tapply(res_GT$est,res_GT$type,function(x)mean(x>obs_est_real["sigma2_I2"]))


save.image("Data/Intermediate/RR_example_results.RData")

