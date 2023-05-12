

rm(list=ls())

run=FALSE


options(mc.cores = parallel::detectCores())
library("beeswarm")

library(squidSim)
library(lme4)
library(parallel)
library(rstan)
library(MCMCglmm)
library(lme4)
rstan_options("auto_write" = TRUE)

wd <- "~/github/null_distributions/"

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
	LMM_stanN <- stan_model(file = paste0(wd,"stan/simple_LMM_normal.stan"))
	LMM_stanU_var <- stan_model(file = paste0(wd,"stan/simple_LMM_uniform_variance.stan"))

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

		stan_dat3 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			cauchy_scale=25)

		stan_mod3 <- sampling(LMM_stan, data=stan_dat3, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_dat4 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			uniform_max=5)

		stan_mod4 <- sampling(LMM_stanU, data=stan_dat4, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_dat5 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			uniform_max=25)

		stan_mod5 <- sampling(LMM_stanU, data=stan_dat5, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_dat6 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			uniform_max=25)

		stan_mod6 <- sampling(LMM_stanU_var, data=stan_dat6, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		stan_dat7 <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"],
			normal_scale=1)

		stan_mod7 <- sampling(LMM_stanN, data=stan_dat7, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)


		prior <- list(G = list(g1=list(V = 1e-16, nu = -2)),R = list(V = 1e-16, nu = -2))

		mcmc_mod <- MCMCglmm(y~1,random=~ID,data=dat[[i]], prior=prior, verbose=FALSE)
		mcmc_post <- mcmc_mod$VCV[,"ID"]
		
		modF <- lmer(y~1+(1|ID),data=dat[[i]])

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
				prior="C25",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod3)[1:4]
			),
			data.frame(
				prior="U5",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod4)[1:4]
			),
			data.frame(
				prior="U25",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod5)[1:4]
			),			
			data.frame(
				prior="U25var",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod6)[1:4]
			),	
			data.frame(
				prior="N1",
				type=c("mode0.1","mode1","median","mean"),
				estimate=stan_out(stan_mod7)[1:4]
			),
			data.frame(
				prior="REML",
				type="freq",
				estimate=as.numeric(summary(modF)$varcor)
			)
			,
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
	out<-prior_sim(n_pop=100, ICC=0.2, N_group=80, N_within=2, mc.cores=8)
	save(out, file=paste0(wd,"Data/Intermediate/prior_impact2.Rdata"))
}
load(paste0(wd,"Data/Intermediate/prior_impact2.Rdata"))

# out <- mclapply(1:500,function(i){
	
# 	modF <- lmer(y~1+(1|ID),data=dat[[i]])

# 	rbind(out[[i]], data.frame(prior="REML",type="freq",estimate=as.numeric(summary(modF)$varcor)))
# }, mc.cores=6)




out2<-do.call(rbind,out)

out2$prior2 <- factor(out2$prior, levels=c("REML", "I", "U25var", "C2", "C5", "C25", "N1", "U5", "U25"))


tail(out2,20)

aggregate(estimate~prior+type,out2,mean)

prior_plot <- function(stat, ...){
	dat <- rbind(subset(out2,type=="freq"),subset(out2,type==stat))
	beeswarm(estimate~prior2,dat, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("REML","Improper","U(0,25)","C(0,2)","C(0,5)","C(0,25)","N(0,1)","U(0,5)","U(0,25)"), ylab="Estimate",... )
	abline(h=0.2, col="red")
	points(aggregate(estimate~prior2,dat,mean)$estimate, cex=1.5, pch=19, col=c("blue",rep("purple",2),rep("orange",6)))
}

setEPS()
pdf(paste0(wd,"Figures/FigSM_prior.pdf"), height=11, width=9)
{
# par(mfrow=c(2,2))
	par(mfrow=c(4,1))
	prior_plot(stat="mode0.1",main="Mode: scale= 0.1")
	prior_plot(stat="mode1",main="Mode: scale= 1")
	prior_plot(stat="mean",main="Mean")
	prior_plot(stat="median",main="Median")


}
dev.off()




plot(subset(out2,type=="mode0.1" & prior=="I")$estimate,subset(out2,type== "freq")$estimate)
plot(subset(out2,type=="mode1" & prior=="I")$estimate,subset(out2,type== "freq")$estimate)


