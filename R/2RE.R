
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
library(MCMCglmm)
data(BTdata)
rstan_options("auto_write" = TRUE)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

LMM_stan <- stan_model(file = paste0(wd,"stan/LMM_normal_2RF.stan"))

data <- BTdata[,c("tarsus","dam", "fosternest")]
names(data) <- c("y","group1","group2")

out_mod<-gaussian_mods2(data)

out_mod$summary
head(out_mod$post)

N_sim=100

sim_dat_both<-simulate_population(
	data_structure= data[,c("group1","group2")],
	parameters= list(residual=list(vcov=median(rowSums(out_mod$post)))),
	n_pop=N_sim
)

sim_mod_both <- mclapply(get_population_data(1:N_sim, list=TRUE),function(j){ 
	gaussian_mods2(sim_dat_both[[j]])$summary
	cat(j," ")
},mc.cores=8)


sim_dat_group1<-simulate_population(
	data_structure= data[,c("group1","group2")],
	parameters= list(residual=list(vcov=median(rowSums(out_mod$post)))),
	n_pop=N_sim
)

sim_dat_group2<-simulate_population(
	data_structure= data[,c("group1","group2")],
	parameters= list(residual=list(vcov=median(rowSums(out_mod$post)))),
	n_pop=N_sim
)


null <- sapply(bootstrap, function(x)x["median","group1"])
	p_func(out_mod$summary["median","group1"],null)
hist(null)

