
rm(list=ls())

devtools::load_all("~/github/squid/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))


set.seed(20220216)
squid_dat <- simulate_population(
	data_structure= make_structure(paste0("ID(200)"),repeat_obs=5),
	parameters= list(
		ID=list(
			names=c("intercepts","slopes"),
			vcov=matrix(c(0.2,0.05,0.05,0.05),2,2),
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
	N_pop=1
)

dat <- get_population_data(squid_dat, list=TRUE)

dat1<-dat[[1]]
var(dat1$y)

plot(y~environment,dat1)
for(i in unique(dat1[,ID])){
	dats <- subset(dat1,ID==i)
	abline(lm(y~environment,dats))
}

mod1 <- lmer(y~environment + (environment|ID),dat1)
mod2 <- lmer(y~environment + (0+environment|ID) + (1|ID),dat1)
mod3 <- lmer(y~environment + (1|ID),dat1)
mod4 <- lm(y~environment,dat1)
anova(mod1,mod2)
anova(mod3,mod2)
anova(mod3,mod4)

summary(mod1)
summary(mod2)



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

		null_mod <- sampling(RR_stan, data=null_dat, chains=1,iter=4000, warmup=1000, pars=c("sigma2_int","sigma2_slope"),, refresh=0)
		stan_out_RR(null_mod)

	})
	cat(i, " ")
	list(param = c(sim_var=0.2, N_group=80, N_within=4), data=dat[[i]],actual=actual,null=null)

},mc.cores=6)

null1 <- sapply(null)

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





set.seed(20211129)
squid_dat<-simulate_population(
	data_structure= make_structure(paste0("ID(80)"),repeat_obs=4),
	parameters= list( ID=list(vcov=0.01, beta=0), residual=list(vcov=0.01, beta=0)),
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

save(results, file=paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))



