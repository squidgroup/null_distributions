
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

set.seed(20220216)
squid_dat <- simulate_population(
	data_structure= make_structure(paste0("ID(80)"),repeat_obs=4),
	parameters= list(
		ID=list(
			names=c("intercepts","slopes"),
			vcorr=matrix(c(0.2,0,0,0.1),2,2),
			beta=c(1,0)
		),
		observation=list(
			names = "environment",
			beta = 0
		), 
		interactions = list(
			names=c("slopes:environment")
		), 
		residual=list(
			vcov=0.7
		)
	),
	N_pop=100
)

dat <- get_population_data(squid_dat, list=TRUE)

# dat1<-dat[[1]]

RR_stan <- stan_model(file = paste0(wd,"stan/rr_LMM.stan"))

results <- mclapply(1:length(dat),function(i){
	# i=1
	# modF <- suppressMessages(lmer(y~environment + (environment|ID),dat[[i]]))
	# modF_null <- glm(y~1,dat[[i]], family="binomial")
	# LRT <- anova(modF,modF_null)$P[2]/2
	# i=1
	stan_dat <- list(
		N = nrow(dat[[i]]),
		N_ID = length(unique(dat[[i]][,ID])),
		y = dat[[i]][,y],
		x_1 = dat[[i]][,environment],
		ID = dat[[i]][,ID])

	stan_mod <- sampling(RR_stan, data=stan_dat, chains=1,iter=4000, warmup=1000, pars=c("sigma2_int","sigma2_slope"), refresh=0)
	# plot(stan_mod)
	# actual<-c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod), LRT=LRT)
	actual<-stan_out_RR(stan_mod)
# traceplot(stan_mod)
	null <- lapply(1:100,function(j){
		perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

		# null_modF <-suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]], REML=FALSE))

		null_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			x_1 = dat[[i]][,environment],
			ID = perm_ID
		)

		null_mod <- sampling(RR_stan, data=null_dat, chains=1,iter=3000, warmup=1000, pars=c("sigma2_int","sigma2_slope"), refresh=0)
		stan_out_RR(null_mod)

		# perm_y <- sample(dat[[i]][,y], replace=FALSE)

		# null_dat <- list(
		# 	N = nrow(dat[[i]]),
		# 	N_ID = length(unique(dat[[i]][,ID])),
		# 	y = perm_y,
		# 	x_1 = dat[[i]][,environment],
		# 	ID = dat[[i]][,ID]
		# )

		# null_mod <- sampling(RR_stan, data=null_dat, chains=1,iter=2000, warmup=1000, pars=c("sigma2_int","sigma2_slope"), refresh=0)
		# stan_out_RR(null_mod)

	})
	cat(i, " ")
	null1 <- do.call(rbind,lapply(null,function(x) x[1,]))
	null2 <- do.call(rbind,lapply(null,function(x) x[2,]))

	list(param = c(sim_var1=0.2,sim_var2=0.1, N_group=80, N_within=4), data=dat[[i]],actual=actual,null1=null1,null2=null2)

},mc.cores=6)


save(results, file=paste0(wd,"Data/Intermediate/RR_sim.Rdata"))




set.seed(20220226)
squid_dat <- simulate_population(
	data_structure= make_structure(paste0("ID(80)"),repeat_obs=4),
	parameters= list(
		ID=list(
			names=c("intercepts","slopes"),
			vcorr=matrix(c(0.2,0,0,0.1),2,2),
			beta=c(0,0)
		),
		observation=list(
			names = "environment",
			beta = 0
		), 
		interactions = list(
			names=c("slopes:environment"),
			beta=c(0,0)
		), 
		residual=list(
			vcov=1
		)
	),
	N_pop=100
)

dat <- get_population_data(squid_dat, list=TRUE)

# dat1<-dat[[1]]

RR_stan <- stan_model(file = paste0(wd,"stan/rr_LMM.stan"))

results <- mclapply(1:length(dat),function(i){

	stan_dat <- list(
		N = nrow(dat[[i]]),
		N_ID = length(unique(dat[[i]][,ID])),
		y = dat[[i]][,y],
		x_1 = dat[[i]][,environment],
		ID = dat[[i]][,ID])

	stan_mod <- sampling(RR_stan, data=stan_dat, chains=1,iter=4000, warmup=1000, pars=c("sigma2_int","sigma2_slope"), refresh=0)

	actual<-stan_out_RR(stan_mod)

	null <- lapply(1:100,function(j){
		perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

		null_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			x_1 = dat[[i]][,environment],
			ID = perm_ID
		)

		null_mod <- sampling(RR_stan, data=null_dat, chains=1,iter=3000, warmup=1000, pars=c("sigma2_int","sigma2_slope"), refresh=0)
		stan_out_RR(null_mod)

	})
	cat(i, " ")
	null1 <- do.call(rbind,lapply(null,function(x) x[1,]))
	null2 <- do.call(rbind,lapply(null,function(x) x[2,]))

	list(param = c(sim_var1=0.2,sim_var2=0.1, N_group=80, N_within=4), data=dat[[i]],actual=actual,null1=null1,null2=null2)

},mc.cores=6)

save(results, file=paste0(wd,"Data/Intermediate/RR_sim_null.Rdata"))
