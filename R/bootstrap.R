
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

n_pop=1000
samples <- as.matrix(expand.grid(ID=c(20,40,80,160),observation=c(2,4)))
ds <- make_structure("ID(500)",repeat_obs=20)

set.seed(20220617)
squid_dat_0 <- simulate_population(
	data_structure= ds,
	parameters= list(residual=list(vcov=1)),
	n_pop=n_pop,
	sample_type="nested",
	sample_param=samples
)

squid_dat_0.2 <- simulate_population(
	data_structure= ds,
	parameters= list( ID=list(vcov=0.2), residual=list(vcov=0.8)),
	n_pop=n_pop,
	sample_type="nested",
	sample_param=samples
)

squid_dat_0.4 <- simulate_population(
	data_structure= ds,
	parameters= list( ID=list(vcov=0.4), residual=list(vcov=0.6)),
	n_pop=n_pop,
	sample_type="nested",
	sample_param=samples
)

save(squid_dat_0,squid_dat_0.2,squid_dat_0.4, file=paste0(wd,"Data/Intermediate/gaussian_sim_data.Rdata"))

dat <- get_sample_data(squid_dat_0, sample_set = 1, list=TRUE)

LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))


dist_M <- mclapply(1:N_pop,function(i){
	},mc.cores=8)
	return(dist_M)

perm_boot <- function(data, N_perm, N_boot){

		modF <- suppressMessages(lmer(y~1+(1|ID),dat[[i]]))

		# i=1
		stan_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID],
			cauchy_scale=2)

		stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID","sigma2_E"), refresh=0)

		actual <- c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod))

		null <- do.call(rbind,lapply(1:N_perm,function(j){
			perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

			null_modF <- suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]]))

			null_dat <- list(
				N = nrow(dat[[i]]),
				N_ID = length(unique(dat[[i]][,ID])),
				y = dat[[i]][,y],
				ID = perm_ID,
				cauchy_scale=2
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

		}))

		all_var <- median(extract(stan_mod)$sigma2_ID + extract(stan_mod)$sigma2_E)

		bs_dat_all<-simulate_population(
			data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
			parameters= list(residual=list(vcov=all_var)),
			n_pop=N_boot
		)
		bs_dat_all2 <- get_population_data(bs_dat_all, list=TRUE)

		bootstrap <- do.call(rbind,lapply(bs_dat_all2,function(j){
			
			bs_modF <- suppressMessages(lmer(y~1+(1|ID),j))

			bs_dat <- list(
				N = nrow(j),
				N_ID = length(unique(j[,ID])),
				y = j[,y],
				ID = j[,ID],
				cauchy_scale=2
			)

			bs_mod <- sampling(LMM_stan, data=bs_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(bs_modF)$varcor),stan_out(bs_mod))

		}))
	

		cat(i, " ")
		list(param = c(ICC=ICC, N_group=N_group, N_within=N_within), data=dat[[i]],actual=actual,null=null, bootstrap=bootstrap)



}




perm_sim <- function(N_pop, N_perm, N_boot, ICC, N_group, N_within,mc.cores=4){

	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=if(ICC==0){list(vcov=0.01, beta=0)}else{list(vcov=ICC)}, residual=list(vcov=1-ICC)),
		n_pop=N_pop
		)

	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))

	dist_M <- mclapply(1:N_pop,function(i){
		modF <- suppressMessages(lmer(y~1+(1|ID),dat[[i]]))

		# i=1
		stan_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID],
			cauchy_scale=2)

		stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID","sigma2_E"), refresh=0)

		actual <- c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod))

		null <- do.call(rbind,lapply(1:N_perm,function(j){
			perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

			null_modF <- suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]]))

			null_dat <- list(
				N = nrow(dat[[i]]),
				N_ID = length(unique(dat[[i]][,ID])),
				y = dat[[i]][,y],
				ID = perm_ID,
				cauchy_scale=2
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

		}))

		all_var <- median(extract(stan_mod)$sigma2_ID + extract(stan_mod)$sigma2_E)

		bs_dat_all<-simulate_population(
			data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
			parameters= list(residual=list(vcov=all_var)),
			n_pop=N_boot
		)
		bs_dat_all2 <- get_population_data(bs_dat_all, list=TRUE)

		bootstrap <- do.call(rbind,lapply(bs_dat_all2,function(j){
			
			bs_modF <- suppressMessages(lmer(y~1+(1|ID),j))

			bs_dat <- list(
				N = nrow(j),
				N_ID = length(unique(j[,ID])),
				y = j[,y],
				ID = j[,ID],
				cauchy_scale=2
			)

			bs_mod <- sampling(LMM_stan, data=bs_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(bs_modF)$varcor),stan_out(bs_mod))

		}))
	

		cat(i, " ")
		list(param = c(ICC=ICC, N_group=N_group, N_within=N_within), data=dat[[i]],actual=actual,null=null, bootstrap=bootstrap)

	},mc.cores=mc.cores)
	return(dist_M)
}

sim_dat <- perm_sim(N_pop=10, N_perm=100,N_boot=100, ICC=0.2, N_group=80, N_within=4,mc.cores=8)

sapply(sim_dat, function(x){
	c(freq = p_func(x$actual["freq"],x$null[,"freq"]),
		freq_b = p_func(x$actual["freq"],x$bootstrap[,"freq"]),
		perm = p_func(x$actual["median"],x$null[,"median"]),
		boot = p_func(x$actual["median"],x$bootstrap[,"median"])
		)

})


#N_pop=10; N_perm=10;N_boot=10; ICC=0.2; N_group=80; N_within=4;mc.cores=8


