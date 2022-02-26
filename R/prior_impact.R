
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

stan_out <- function(model){
	out <- c(MCMCglmm::posterior.mode(coda::as.mcmc(as.data.frame(extract(model)))[,"sigma2_ID"]), 
		median(extract(model)$sigma2_ID),
		summary(model)$summary["sigma2_ID",1])
	names(out) <- c("mode","median","mean")
	out	
}

prior_sim <- function(N_pop, ICC, N_group, N_within, mc.cores=4){
#N_pop=100; ICC=0.2; N_group=80; N_within=2; mc.cores=8
	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=if(ICC==0){list(vcov=0.01, beta=0)}else{list(vcov=ICC)}, residual=list(vcov=1-ICC)),
		N_pop=N_pop
		)
	
	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))
	LMM_stanU <- stan_model(file = paste0(wd,"stan/simple_LMM_uniform.stan"))

#i=1
	dist_M <- mclapply(1:N_pop,function(i){
		modF <- suppressMessages(lmer(y~1+(1|ID),dat[[i]], REML=TRUE))
		
		stan_dat1 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID],
			cauchy_scale=2)

		stan_mod1 <- sampling(LMM_stan, data=stan_dat1, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
		
		stan_dat2 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID],
			cauchy_scale=5)

		stan_mod2 <- sampling(LMM_stan, data=stan_dat2, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_mod3 <- sampling(LMM_stanU, data=stan_dat1, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		cat(i, " ")

		rbind(
			data.frame(
				prior=NA,
				type="freq",
				estimate=as.numeric(summary(modF)$varcor)
			),
			data.frame(
				prior="C2",
				type=c("mode","median","mean"),
				estimate=stan_out(stan_mod1)
			),
			data.frame(
				prior="C5",
				type=c("mode","median","mean"),
				estimate=stan_out(stan_mod2)
			),
			data.frame(
				prior="U",
				type=c("mode","median","mean"),
				estimate=stan_out(stan_mod3)
			)
		)

	
	},mc.cores=mc.cores)
	return(dist_M)
}


out<-prior_sim(N_pop=100, ICC=0.2, N_group=80, N_within=2, mc.cores=8)

out2<-do.call(rbind,out)

boxplot(mode~prior,out2)

library("beeswarm")
	beeswarm(estimate~prior+type,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")


par(mfrow=c(2,2))
	beeswarm(mode~prior,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")
		beeswarm(freq~prior,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")
			beeswarm(mean~prior,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")
				beeswarm(median~prior,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")



lapply(out,funtion(x) x)
