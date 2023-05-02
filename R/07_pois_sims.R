
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
ds <- make_structure(paste0("ID(80)"),repeat_obs=2)

# between group variances 
ICCs <- c(0, 0.1, 0.2, 0.4)
intercept <- 2
total_var <- 0.2
## match median mean and variances from Pick et al 2023

# (0.449*8.41)^2
#  exp2lat(8.41,14)
# lat2exp(2,0.2)

# n=100000
# x<-rpois(n,rlnorm(n,1.5,sqrt(0.1)))
# hist(x)
# mean(x);var(x)

# import stan code for mixed model
pois_stan <- stan_model(file = paste0(wd,"stan/simple_pois_GLMM.stan"))


set.seed("20230427")
# ICC<-0
# i=1
for(ICC in ICCs){
	cat(ICC,"\n")
	squid_dat<-simulate_population(
		data_structure= ds,
		parameters= list( 
			intercept= intercept,
			ID=list(vcov=ICC*total_var), 
			residual=list(vcov=total_var*(1-ICC))
		),
		family="poisson", link="log",
		n_pop=n_pop
	)

	dat <- get_population_data(squid_dat,list=TRUE)
	
	out <- mclapply(1:length(dat),function(i){
		cat(i, " ")
		out_mod<-pois_mods(dat[[i]])
		out_mod[["param"]]<- c(pop=i,ICC=ICC, N_group=80, N_within=2)
		out_mod
	},mc.cores=4)
	save(out, file=paste0(wd,"Data/Intermediate/pois_sims_",ICC,".Rdata"))
		cat("\n\n")
}

	
