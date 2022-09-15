
rm(list=ls())

options(mc.cores = parallel::detectCores())

# library(lme4)
library(squidSim)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

GLMM_stan <- stan_model(file = paste0(wd,"stan/simple_GLMM.stan"))

set.seed(20211128)

for(ICC in c(0,0.2)){
	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(80)"),repeat_obs=4),
		parameters= list( 
			ID=list(vcov=ICC), 
			residual=list(vcov=0)
		),
		family="binomial", link="logit",
		n_pop=100
	)

	dat <- get_population_data(squid_dat, list=TRUE)

	results <- mclapply(1:length(dat),function(i){
		# i=1

		stan_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,"ID"])),
			y = dat[[i]][,"y"],
			ID = dat[[i]][,"ID"])

		stan_mod <- sampling(GLMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		actual<-stan_out(stan_mod)

	## permutation method
		perm_out <- do.call(rbind,lapply(1:100,function(j){

			perm_dat <- list(
				N = nrow(dat[[i]]),
				N_ID = length(unique(dat[[i]][,"ID"])),
				y = dat[[i]][,"y"],
				ID = sample(dat[[i]][,"ID"], replace=FALSE)
			)

			perm_mod <- sampling(GLMM_stan, data=perm_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

			stan_out(perm_mod)

		}))

	## bootstrap method
		bs_dat_all <- simulate_population(
			data_structure = squid_dat$data_structure,
			parameters = list(
				ID=list(vcov=0),
				residual=list(vcov=0)
			),
			family="binomial", link="logit",
			n_pop=100
		)
		bs_dat <- get_population_data(bs_dat_all, list=TRUE)

		boot_out <- do.call(rbind,lapply(bs_dat,function(x){
			
			boot_dat <- list(
				N = nrow(x),
				N_ID = length(unique(x[,"ID"])),
				y = x[,"y"],
				ID = x[,"ID"]
			)
			boot_mod <- sampling(GLMM_stan, data=boot_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

			stan_out(boot_mod)
		}))

		cat(i, " ")
		list(
			param = c(sim_var=ICC, N_group=80, N_within=4),
			actual=actual,
			perm=perm_out,
			boot=boot_out
		)
	},mc.cores=8)
	assign(paste0("results",ICC),results)
	cat("\n\n")
}
save(results0,results0.2, file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))










# actual <- as.data.frame(
# 	t(sapply(results, function(x) c(x$param,x$actual[1:6])))
# )

# long<-Stack(actual,c("mean","median","mode","freq"), group.name="type",value.name = "estimate")

# 	boxplot(estimate~N_within+type,long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)
# 	abline(h=0.2,col="red",lwd=2)

# colMeans(actual[,c("mean","median","mode","freq")])
# apply(actual[,c("mean","median","mode","freq")],2,median)




# p_perm <- as.data.frame(
# 	t(sapply(results, function(x) c(x$param,
# 		freq = p_func(x$actual["freq"],x$null[,"freq"]),
# 		mean=p_func(x$actual["mean"],x$null[,"mean"]),
# 		median=p_func(x$actual["median"],x$null[,"median"]),
# 		mode= p_func(x$actual["mode"],x$null[,"mode"]),
# 		x$actual["LRT"]
# 		)))
# )

# p_long<-Stack(p_perm,c("mean","median","mode","freq","LRT"), group.name="type",value.name = "estimate")

# 	boxplot(estimate~N_within+type,p_long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)


# apply(p_perm[,c("mean","mode","freq","LRT")],2, function(x) mean(x<0.05))

# boxplot(sapply(results, function(x) c(x$actual["LRT"])))



