
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

## number of populations per parameter set
n_pop=1000

## population data_structure
ds <- make_structure("ID(500)",repeat_obs=20)

## between and within group sample sizes
samples <- as.matrix(expand.grid(ID=c(20,40,80,160),observation=c(2,4)))

LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))


######
## No between group variance
######

set.seed("202206171")

ICC=0

squid_dat_0 <- simulate_population(
	data_structure= ds,
	parameters= list(residual=list(vcov=1-ICC)),
	n_pop=n_pop,
	sample_type="nested",
	sample_param=samples
)

for (j in 1:nrow(squid_dat_0$sample_param)){
	dat <- get_sample_data(squid_dat_0, sample_set = j, list=TRUE)
	param <- squid_dat_0$sample_param[j,]
		cat("\n set:",param, "\n")
	out <- mclapply(1:length(dat),function(i){
		cat(i, " ")
		c(gaussian_mods(dat[[i]]),ICC=ICC, N_group=param[1], N_within=param[2])
	},mc.cores=8)
	save(out, file=paste0(wd,"Data/Intermediate/gaus_sims_",ICC,"_",param[1],"_",param[2],".Rdata"))

}



######
## Between group variance
######

set.seed("202206172")

for(ICC in c(0.2,0.4)){

	squid_dat_0 <- simulate_population(
		data_structure= ds,
		parameters= list( ID=list(vcov=ICC), residual=list(vcov=1-ICC)),
		n_pop=n_pop,
		sample_type="nested",
		sample_param=samples
	)

	for (j in 1:nrow(squid_dat_0$sample_param)){
		dat <- get_sample_data(squid_dat_0, sample_set = j, list=TRUE)
		param <- squid_dat_0$sample_param[j,]
			cat("\n set:",param, "\n")
		out <- mclapply(1:length(dat),function(i){
			cat(i, " ")
			c(gaussian_mods(dat[[i]]),ICC=ICC, N_group=param[1], N_within=param[2])
		},mc.cores=8)
		save(out, file=paste0(wd,"Data/Intermediate/gaus_sims_",ICC,"_",param[1],"_",param[2],".Rdata"))
	}

}

	
