
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

## number of populations per parameter set
n_pop=500

## population data_structure - make large enough, so the samples are unlikely to contain much overlap
ds <- make_structure("ID(300)",repeat_obs=10)

## between and within group sample sizes
samples <- as.matrix(expand.grid(ID=c(20,40,80),observation=c(2,4)))

# between group variances 
ICCs <- c(0,0.1,0.2,0.4)

# import stan code for mixed model
LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))


set.seed("20220617")

for(ICC in ICCs){

	squid_dat <- simulate_population(
		data_structure= ds,
		parameters= list( ID=list(beta=sqrt(ICC)), residual=list(vcov=1-ICC)),
		n_pop=n_pop,
		sample_type="nested",
		sample_param=samples, verbose=TRUE
	)

	for (j in 1:nrow(squid_dat$sample_param)){
		dat <- get_sample_data(squid_dat, sample_set = j, list=TRUE)
		param <- squid_dat$sample_param[j,]
			cat("\n set",j,":",ICC,param, "\n")
		out <- mclapply(1:length(dat),function(i){
			cat(i, " ")
			out_mod<-gaussian_mods(dat[[i]])
			out_mod[["param"]]<- c(pop=i,ICC=ICC, N_group=param[1], N_within=param[2])
			out_mod
		},mc.cores=8)
		save(out, file=paste0(wd,"Data/Intermediate/gaus_sims_",ICC,"_",param[1],"_",param[2],".Rdata"))
	}

}

	
