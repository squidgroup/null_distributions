
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
sim_dat_both_2 <- get_population_data(sim_dat_both, list=TRUE)

sim_mod_both <- mclapply(1:N_sim,function(j){ 
	cat(j," ")
	gaussian_mods2(sim_dat_both_2[[j]])$summary
},mc.cores=8)


sim_dat_group1<-simulate_population(
	data_structure= data[,c("group1","group2")],
	parameters= list(
		group2=list(vcov=out_mod$summary["median","group2"]),
		residual=list(vcov=median(rowSums(out_mod$post[,c(1,3)])))
	),
	n_pop=N_sim
)
sim_dat_group1_2<- get_population_data(sim_dat_group1, list=TRUE)

sim_mod_group1 <- mclapply(1:N_sim,function(j){ 
	cat(j," ")
	gaussian_mods2(sim_dat_group1_2[[j]])$summary
},mc.cores=8)


sim_dat_group2<-simulate_population(
	data_structure= data[,c("group1","group2")],
	parameters= list(
		group1=list(vcov=out_mod$summary["median","group1"]),
		residual=list(vcov=median(rowSums(out_mod$post[,c(2,3)])))
	),
	n_pop=N_sim
)
sim_dat_group2_2<-get_population_data(sim_dat_group2, list=TRUE)

sim_mod_group2 <- mclapply(1:N_sim,function(j){ 
	cat(j," ")
	gaussian_mods2(sim_dat_group2_2[[j]])$summary
},mc.cores=8)


null_both1 <- sapply(sim_mod_both, function(x)x["median","group1"])
null_both2 <- sapply(sim_mod_both, function(x)x["median","group2"])

null_1 <- sapply(sim_mod_group1, function(x)x["median","group1"])
null_2 <- sapply(sim_mod_group2, function(x)x["median","group2"])

p_func(out_mod$summary["median","group1"],null_both1)
p_func(out_mod$summary["median","group2"],null_both2)
p_func(out_mod$summary["median","group1"],null_1)
p_func(out_mod$summary["median","group2"],null_2)

nulls <- data.frame(null=c(null_both1,null_both2,null_1,null_2),dist=rep(c("1both","2both","1","2"),each=100))

beeswarm::beeswarm(null~dist,nulls, pch=19, cex=0.5, col=alpha(1,0.2),method = "compactswarm",corral="wrap")

