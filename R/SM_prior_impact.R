
rm(list=ls())

run=FALSE


options(mc.cores = parallel::detectCores())
library("beeswarm")

library(squidSim)
library(lme4)
library(parallel)
library(rstan)
library(MCMCglmm)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

prior_sim <- function(n_pop, ICC, N_group, N_within, mc.cores=4){
#n_pop=100; ICC=0.2; N_group=80; N_within=2; mc.cores=8
	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=list(vcov=ICC), residual=list(vcov=1-ICC)),
		n_pop=n_pop
		)
	
	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))
	LMM_stanU <- stan_model(file = paste0(wd,"stan/simple_LMM_uniform.stan"))

#i=1
	dist_M <- mclapply(1:n_pop,function(i){
		
		stan_dat1 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			cauchy_scale=2)

		stan_mod1 <- sampling(LMM_stan, data=stan_dat1, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
		
		stan_dat2 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			cauchy_scale=5)

		stan_mod2 <- sampling(LMM_stan, data=stan_dat2, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_mod3 <- sampling(LMM_stanU, data=stan_dat1, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		prior <- list(G = list(g1=list(V = 1e-16, nu = -2)),R = list(V = 1e-16, nu = -2))

		mcmc_mod <- MCMCglmm(y~1,random=~ID,data=dat[[i]], prior=prior, verbose=FALSE)
		mcmc_post <- mcmc_mod$VCV[,"ID"]
		
		## maybe add in freq

		cat(i, " ")

		rbind(
			data.frame(
				prior="C2",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod1)[1:4]
			),
			data.frame(
				prior="C5",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod2)[1:4]
			),
			data.frame(
				prior="U",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod3)[1:4]
			),
			data.frame(
				prior="I",
				type=c("mode0.1","mode1","median","mean"),
				estimate=c(post_mode(mcmc_post,adjust=0.1), post_mode(mcmc_post,adjust=1),  median(mcmc_post), mean(mcmc_post))
			)
		)


	
	},mc.cores=mc.cores)
	return(dist_M)
}

if(run){
	set.seed(20221005)
	out<-prior_sim(n_pop=500, ICC=0.2, N_group=80, N_within=2, mc.cores=8)
	save(out, file=paste0(wd,"Data/Intermediate/prior_impact.Rdata"))
}
load(paste0(wd,"Data/Intermediate/prior_impact.Rdata"))


out2<-do.call(rbind,out)


prior_plot <- function(stat, ...){
	beeswarm(estimate~prior,subset(out2,type==stat), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Improper","Uniform(0,2)"), ylab="Estimate",... )
	abline(h=0.2, col="red")
	points(aggregate(estimate~prior,subset(out2,type==stat),mean)$estimate, cex=1.5, pch=19, col="orange")
}

setEPS()
pdf(paste0(wd,"Figures/FigSM_prior.pdf"), height=9, width=9)
{
par(mfrow=c(2,2))
	prior_plot(stat="mode0.1",main="Mode: scale= 0.1")
	prior_plot(stat="mode1",main="Mode: scale= 1")
	prior_plot(stat="mean",main="Mean")
	prior_plot(stat="median",main="Median")


}
dev.off()

	beeswarm(estimate~prior+type,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")

setEPS()
pdf(paste0(wd,"Figures/FigSM_prior.pdf"), height=9, width=9)
{
par(mfrow=c(2,2))
	
	# ICC=0.2, N_group=80, N_within=2
	beeswarm(estimate~prior,subset(out2,type=="mode1"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Mode: scale= 1")
	abline(h=0.2, col="red")
	points(aggregate(estimate~prior,subset(out2,type=="mode1"),mean)$estimate, cex=1.5, pch=19, col="orange")

	beeswarm(estimate~prior,subset(out2,type=="mean"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Mean")
	abline(h=0.2, col="red")
	points(aggregate(estimate~prior,subset(out2,type=="mean"),mean)$estimate, cex=1.5, pch=19, col="orange")

	beeswarm(estimate~prior,subset(out2,type=="median"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Median")
	abline(h=0.2, col="red")
	points(aggregate(estimate~prior,subset(out2,type=="median"),mean)$estimate, cex=1.5, pch=19, col="orange")

}
dev.off()


