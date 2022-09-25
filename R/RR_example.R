library(squidSim)
library(arm)
library(rstan)
library(parallel)
rm(list=ls())
setwd("/home/yi/Desktop/Perm_RRExample")
source("R/RR_functions.R")

#Simulated example
##Generate "real" simulated data setlibrary(squidSim)
library(arm)
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
obs_md_sim <- stan("Stan/RR_example.stan", data = stan_data, pars = params, chains = 3, iter = 5500,  warmup = 500, thin = 5, cores=3)
obs_est_sim<-round(summary(obs_md_sim)$summary[,1],3)

##Permutation tests
n_perm<-100
t1<-system.time({  
l_rand_ID_data<-ran_ID(data, n_perm)
res_ID<-unlist(mclapply(l_rand_ID_data, mod, mc.cores=8))
})

l_rand_y_data<-ran_y(data, n_perm)
res_y<-unlist(mclapply(l_rand_y_data, mod, mc.cores=8))

l_rand_x_data<-ran_x(data, n_perm)
res_x<-unlist(mclapply(l_rand_x_data, mod, mc.cores=8))

l_rand_x_ID_data<-ran_x_ID(data, n_perm)
res_x_ID<-unlist(mclapply(l_rand_x_ID_data, mod, mc.cores=8))

##Bootstrapping
t2<-system.time({  
squid_data_null <- simulate_population(data_structure=make_structure("individual(300)",repeat_obs=4), n_pop = n_perm,
                                       parameters = list(
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
                                         )),  known_predictors = list(
                                           predictors = data.frame(environment=data[,"environment"]), 
                                           beta = obs_est_sim[2])
)
data_null <- get_population_data(squid_data_null, list=TRUE)
res_boot<-unlist(mclapply(squid_data_null, mod, mc.cores=8))
})

##Results data frame
res<-data.frame(est=c(res_y, res_ID, res_x, res_x_ID, res_boot), 
                type=as.factor(rep(c("perm_y", "perm_ID", "perm_x", "perm_x_ID","boot"), each=n_perm)))


##P-values
table(res_y<obs_est_sim[4])[1]/sum(table(res_y>obs_est_sim[4]))
table(res_ID<obs_est_sim[4])[1]/sum(table(res_ID>obs_est_sim[4]))
table(res_x<obs_est_sim[4])[1]/sum(table(res_x>obs_est_sim[4]))
table(res_x_ID<obs_est_sim[4])[1]/sum(table(res_x_ID>obs_est_sim[4]))
table(res_boot<obs_est_sim[4])[1]/sum(table(res_boot>obs_est_sim[4]))


##Real example
data_GT<-read.csv("Data/data_GT.csv")
stan_data_GT = list(n = nrow(data_GT),
                 n_ind=max(data_GT$individual),
                 individual=data_GT$individual,
                 x=data_GT$environment,
                 y=data_GT$y)

params<-c("beta_0", "beta_1", "sigma2_I1", "sigma2_I2", "sigma2_e")
obs_md_GT <- stan("Stan/RR_example.stan", data = stan_data_GT, pars = params, chains = 3, iter = 5500,  warmup = 500, thin = 30, cores=3)
obs_est_real<-round(summary(obs_md_GT)$summary[,1],3)

##Permutation tests
t3<-system.time({  
l_rand_ID_data_GT<-ran_ID(data_GT, n_perm)
res_ID_GT<-unlist(mclapply(l_rand_ID_data_GT, mod, mc.cores=8))
})

l_rand_y_data_GT<-ran_y(data_GT, n_perm)
res_y_GT<-unlist(mclapply(l_rand_y_data_GT, mod, mc.cores=8))

l_rand_x_data_GT<-ran_x(data_GT, n_perm)
res_x_GT<-unlist(mclapply(l_rand_x_data_GT, mod, mc.cores=8))

l_rand_x_ID_data_GT<-ran_x_ID(data_GT, n_perm)
res_x_ID_GT<-unlist(mclapply(l_rand_x_ID_data_GT, mod, mc.cores=8))

##Bootsraping
t4<-system.time({  
  
squid_data_null_GT <- simulate_population(data_structure=data.frame(individual=data_GT$individual), n_pop = n_perm,
                                       parameters = list(
                                         individual = list(
                                           names = c("ind_int","ind_slope"), 
                                           beta = c(1,0),
                                           vcov = c(obs_est_real[3],0)
                                         ),
                                         residual = list(
                                           vcov = c(obs_est_real[4] + obs_est_real[5])
                                         ),
                                         interactions = list(
                                           names = c("ind_slope:environment"),
                                           beta = c(1)
                                         )),  known_predictors = list(
                                           predictors = data.frame(environment=data_GT[,"environment"]), 
                                           beta = obs_est_real[2])
)


data_null_GT <- get_population_data(squid_data_null, list=TRUE)
res_boot_GT<-unlist(mclapply(data_null_GT, mod, mc.cores=8))
})

res_GT<-data.frame(est=c(res_y_GT, res_ID_GT, res_x_GT, res_x_ID_GT, res_boot_GT), 
                type=as.factor(rep(c("perm_y", "perm_ID", "perm_x", "perm_x_ID","boot"), each=n_perm)))


table(res_y_GT<obs_est_real[4])[1]/sum(table(res_y_GT>obs_est_real[4]))
table(res_ID_GT<obs_est_real[4])[1]/sum(table(res_ID_GT>obs_est_real[4]))
table(res_x_GT<obs_est_real[4])[1]/sum(table(res_x_GT>obs_est_real[4]))
table(res_x_ID_GT<obs_est_real[4])[1]/sum(table(res_x_ID_GT>obs_est_real[4]))
table(res_boot_GT<obs_est_real[4])[1]/sum(table(res_boot_GT>obs_est_real[4]))

save.image("Data/RR_example_results.RData")

