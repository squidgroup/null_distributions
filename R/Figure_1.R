
rm(list=ls())

devtools::load_all("~/github/squid/R")

library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/github/bayes_perm/"

stan_out <- function(model){
	out <- c(MCMCglmm::posterior.mode(as.mcmc(as.data.frame(extract(model)))[,"sigma2_ID"]), summary(model)$summary["sigma2_ID",c(1,4,8,9)])
	names(out) <- c("mode","mean","LCI","UCI","ESS")
	out	
}
p_func <- function(actual,null) mean(actual<null)

set.seed(1255)
squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(80)"),repeat_obs=2),
		parameters= list( ID=list(vcov=0.2), residual=list(vcov=0.8)),
		N_pop=1
		)
dat <- get_population_data(squid_dat)

LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))

stan_dat <- list(
			N = nrow(dat),
			N_ID = length(unique(dat[,ID])),
			y = dat[,y],
			ID = dat[,ID])

stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
summary(stan_mod)
hist(extract(stan_mod)$sigma2_ID, breaks=100)

		actual<-stan_out(stan_mod)

	null <- do.call(rbind,mclapply(1:100,function(j){
			perm_ID <- sample(dat[,ID], replace=FALSE)

			null_modF <-suppressMessages(lmer(y~1+(1|perm_ID),dat, REML=FALSE))

			null_dat <- list(
				N = nrow(dat),
				N_ID = length(unique(dat[,ID])),
				y = dat[,y],
				ID = perm_ID
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

		}))

	hist(null[,"mean"])
p_func(actual["mean"],null[,"mean"])






# prior <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=0.002)))
prior_PE <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alphaV=1000)))
# prior_flat <-list(R = list(V = 1, nu = 0), G=list(G1 = list(V = 1, nu = 0)))
modM<-MCMCglmm(y~1, random=~ID,data=dat)
# modM1<-MCMCglmm(y~1,data=dat)
summary(modM)$Gcov[1]
plot(modM)

nullmodelsVCV_M <- sapply(1:100,function(i){
	dat2<-dat
	dat2$nully <- sample(dat2$y, nrow(dat2))
    nullmodel <- MCMCglmm(nully~1, random=~ID,data=dat2)
     var_outM(nullmodel)
})

nullmodelsVCV_M
hist(nullmodelsVCV_M, breaks=100)
abline(v=summary(modM)$Gcov[1], col="red")

