
rm(list=ls())

run=FALSE

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

prior_sim <- function(N_pop, ICC, N_group, N_within, mc.cores=4){
#N_pop=100; ICC=0.2; N_group=80; N_within=2; mc.cores=8
	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=if(ICC==0){list(vcov=0.01, beta=0)}else{list(vcov=ICC)}, residual=list(vcov=1-ICC)),
		n_pop=N_pop
		)
	
	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))
	LMM_stanU <- stan_model(file = paste0(wd,"stan/simple_LMM_uniform.stan"))

#i=1
	dist_M <- mclapply(1:N_pop,function(i){
		
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
			)
		)

	
	},mc.cores=mc.cores)
	return(dist_M)
}

if(run){
out<-prior_sim(N_pop=500, ICC=0.2, N_group=80, N_within=2, mc.cores=8)

subset(out2)


diffs<-sapply(out,function(x){
 median <- subset(x,type=="median")
 mean <- subset(x,type=="mean")
 mode1 <- subset(x,type=="mode1")
 mode0.1 <- subset(x,type=="mode0.1")

 c(medianC5 = median[median$prior=="C2","estimate"]- median[median$prior=="C5","estimate"],
  medianU = median[median$prior=="C2","estimate"]- median[median$prior=="U","estimate"],
 
  meanC5 = mean[mean$prior=="C2","estimate"]- mean[mean$prior=="C5","estimate"],
  meanU = mean[mean$prior=="C2","estimate"]- mean[mean$prior=="U","estimate"],

	mode1C5 = mode1[mode1$prior=="C2","estimate"]- mode1[mode1$prior=="C5","estimate"],
  mode1U = mode1[mode1$prior=="C2","estimate"]- mode1[mode1$prior=="U","estimate"],
  
  mode0.1C5 = mode0.1[mode0.1$prior=="C2","estimate"]- mode0.1[mode0.1$prior=="C5","estimate"],
  mode0.1U = mode0.1[mode0.1$prior=="C2","estimate"]- mode0.1[mode0.1$prior=="U","estimate"],
  mode0.1C5U = mode0.1[mode0.1$prior=="C5","estimate"]- mode0.1[mode0.1$prior=="U","estimate"]
  )
})

boxplot(t(diffs))
rowMeans(diffs)
apply(diffs,1,var)

aggregate(estimate~prior+type,out2,mean)

save(out, file=paste0(wd,"Data/Intermediate/prior_impact.Rdata"))
}
load(paste0(wd,"Data/Intermediate/prior_impact.Rdata"))

out2<-do.call(rbind,out)
head(out2,20)
library("beeswarm")
	beeswarm(estimate~prior+type,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")

setEPS()
pdf(paste0(wd,"Figures/FigSM_prior.pdf"), height=9, width=9)
{
par(mfrow=c(2,2))
	beeswarm(estimate~prior,subset(out2,type=="mode0.1"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Mode: scale= 0.1")
	beeswarm(estimate~prior,subset(out2,type=="mode1"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Mode: scale= 1")
	beeswarm(estimate~prior,subset(out2,type=="mean"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Mean")
	beeswarm(estimate~prior,subset(out2,type=="median"), pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Prior", labels=c("Cauchy(0,2)","Cauchy(0,5)","Uniform(0,2)"), ylab="Estimate", main="Median")
}
dev.off()


