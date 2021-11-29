
rm(list=ls())

devtools::load_all("~/github/squid/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))


set.seed(20211128)
squid_dat<-simulate_population(
	data_structure= make_structure(paste0("ID(80)"),repeat_obs=4),
	parameters= list( ID=list(vcov=0.2), residual=list(vcov=0.01, beta=0)),
	family="binomial", link="logit",
	N_pop=100
	)

dat <- get_population_data(squid_dat, list=TRUE)

GLMM_stan <- stan_model(file = paste0(wd,"stan/simple_GLMM.stan"))

results <- mclapply(1:length(dat),function(i){
	# i=1
	modF <- suppressMessages(glmer(y~1+(1|ID),dat[[i]], family="binomial"))
	modF_null <- glm(y~1,dat[[i]], family="binomial")
	LRT <- anova(modF,modF_null)$P[2]/2
	# i=1
	stan_dat <- list(
		N = nrow(dat[[i]]),
		N_ID = length(unique(dat[[i]][,ID])),
		y = dat[[i]][,y],
		ID = dat[[i]][,ID])

	stan_mod <- sampling(GLMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

	actual<-c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod), LRT=LRT)

	null <- do.call(rbind,lapply(1:100,function(j){
		perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

		null_modF <-suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]], REML=FALSE))

		null_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = perm_ID
		)

		null_mod <- sampling(GLMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
		c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

	}))
	cat(i, " ")
	list(param = c(sim_var=0.2, N_group=80, N_within=4), data=dat[[i]],actual=actual,null=null)

},mc.cores=6)

save(results, file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))


actual <- as.data.frame(
	t(sapply(results, function(x) c(x$param,x$actual[1:6])))
)

long<-Stack(actual,c("mean","median","mode","freq"), group.name="type",value.name = "estimate")

	boxplot(estimate~N_within+type,long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)
	abline(h=0.2,col="red",lwd=2)

colMeans(actual[,c("mean","median","mode","freq")])
apply(actual[,c("mean","median","mode","freq")],2,median)




p_perm <- as.data.frame(
	t(sapply(results, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		median=p_func(x$actual["median"],x$null[,"median"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"]),
		x$actual["LRT"]
		)))
)

p_long<-Stack(p_perm,c("mean","median","mode","freq","LRT"), group.name="type",value.name = "estimate")

	boxplot(estimate~N_within+type,p_long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)


apply(p_perm[,c("mean","mode","freq","LRT")],2, function(x) mean(x<0.05))

boxplot(sapply(results, function(x) c(x$actual["LRT"])))
