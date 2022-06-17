
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

n_pop=200

LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))



perm_boot <- function(output, N_perm, N_boot){

	data <- output$data
  ICC=output$ICC
  N_group=output$N_group
  N_within=output$N_within

	perm <- do.call(rbind,lapply(1:N_perm,function(j){
		data[,"ID"] <- sample(data[,"ID"], replace=FALSE)

		c(gaussian_mods(data)$summary)
	}))

	bs_dat_all<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list(residual=list(vcov=median(rowSums(data$post)))),
		n_pop=N_boot
	)

	bootstrap <- do.call(rbind,lapply(get_population_data(bs_dat_all, list=TRUE),function(j){
		
		c(gaussian_mods(j)$summary)

	}))

	list(param = c(ICC=ICC, N_group=N_group, N_within=N_within),output$summary,actual=perm=perm, bootstrap=bootstrap)

}
#	cat(i, " ")




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
