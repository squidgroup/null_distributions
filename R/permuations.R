
rm(list=ls())

devtools::load_all("~/github/squid/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))


perm_sim <- function(N_pop, N_perm, ICC, N_group, N_within,mc.cores=4){

	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=if(ICC==0){list(vcov=0.01, beta=0)}else{list(vcov=ICC)}, residual=list(vcov=1-ICC)),
		N_pop=N_pop
		)

	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))

	dist_M <- mclapply(1:N_pop,function(i){
		modF <- suppressMessages(lmer(y~1+(1|ID),dat[[i]], REML=FALSE))
		modF_null <- lm(y~1,dat[[i]])
		LRT <- anova(modF,modF_null)$P[2]/2
		# i=1
		stan_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID])

		stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		actual <- c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod), LRT=LRT)

		null <- do.call(rbind,lapply(1:N_perm,function(j){
			perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

			null_modF <- suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]], REML=FALSE))
			null_LRT <- anova(null_modF,modF_null)$P[2]/2

			null_dat <- list(
				N = nrow(dat[[i]]),
				N_ID = length(unique(dat[[i]][,ID])),
				y = dat[[i]][,y],
				ID = perm_ID
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod), LRT=null_LRT)

		}))
		cat(i, " ")
		list(param = c(ICC=ICC, N_group=N_group, N_within=N_within), data=dat[[i]],actual=actual,null=null)

	},mc.cores=mc.cores)
	return(dist_M)
}


set.seed(20211129)
system.time({
for(N_within in c(2,4,8)){
	for(ICC in c(0,0.2,0.4)){
		for(N_group in c(20,40,80,160)){		
			cat("\n",ICC,N_group,N_within,"\n")
			sim_dat <- perm_sim(N_pop=100, N_perm=100, ICC=ICC, N_group=N_group, N_within=N_within,mc.cores=6)
			save(sim_dat, file=paste0(wd,"Data/Intermediate/sims_",ICC,"_",N_group,"_",N_within,".Rdata"))
		}
	}
}
})


## prior used
x <- seq(0,1,length.out=10000)
plot(x,dcauchy(x,0,2), ylim=c(0,max(dcauchy(x,0,2))))
