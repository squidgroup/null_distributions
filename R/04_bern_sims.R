
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

## population data_structure - will be the same for all in these sims
ds <- make_structure(paste0("ID(80)"),repeat_obs=4)


# between group variances 
ICCs <- c(0, 0.2, 0.4, 0.8)

# import stan code for mixed model
bern_stan <- stan_model(file = paste0(wd,"stan/simple_GLMM.stan"))


set.seed("20211128")

for(ICC in ICCs){
	cat(ICC,"\n")
	squid_dat<-simulate_population(
		data_structure= ds,
		parameters= list( 
			ID=list(vcov=ICC), 
			residual=list(vcov=0)
		),
		family="binomial", link="logit",
		n_pop=n_pop
	)

	dat <- get_population_data(squid_dat,list=TRUE)
	
	out <- mclapply(1:length(dat),function(i){
		cat(i, " ")
		out_mod<-bern_mods(dat[[i]])
		out_mod[["param"]]<- c(pop=i,ICC=ICC, N_group=80, N_within=4)
		out_mod
	},mc.cores=8)
	save(out, file=paste0(wd,"Data/Intermediate/bern_sims_",ICC,".Rdata"))
		cat("\n\n")
}

	
