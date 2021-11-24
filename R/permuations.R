
rm(list=ls())

devtools::load_all("~/github/squid/R")

options(mc.cores = parallel::detectCores())

library(lme4)

# library(MCMCglmm)
library(parallel)
# library(MCMCpack)
library(rstan)
library(posterior)
library(coda)

stan_out <- function(model){
	out <- c(MCMCglmm::posterior.mode(as.mcmc(as.data.frame(extract(model)))[,"sigma2_ID"]), summary(model)$summary["sigma2_ID",c(1,4,8,9)])
	names(out) <- c("mean","mode","LCI","UCI","ES")
	out	
}

wd <- "~/github/bayes_perm"

# N_pop <- 100
# N_perm <- 100
# ICC <- c(0,0.2,0.4,0.8)
# N_cluster <- c(20,40,80,160)
# N_within <- c(2,4,8,16)

perm_sim <- function(N_pop, N_perm, ICC, N_group, N_within,mc.cores=4){
	squid_dat<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list( ID=list(vcov=ICC), residual=list(vcov=1-ICC)),
		N_pop=N_pop
		)

	dat <- get_population_data(squid_dat, list=TRUE)

	LMM_stan <- stan_model(file = paste0(wd,"/stan/simple_LMM.stan"))

	dist_M <- mclapply(1:N_pop,function(i){
		modF <-suppressMessages(lmer(y~1+(1|ID),dat[[i]], REML=FALSE))
		# i=1
		stan_dat <- list(
			N = nrow(dat[[i]]),
			N_ID = length(unique(dat[[i]][,ID])),
			y = dat[[i]][,y],
			ID = dat[[i]][,ID])

		stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		actual<-c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod))

		null <- do.call(rbind,lapply(1:N_perm,function(j){
			perm_ID <- sample(dat[[i]][,ID], replace=FALSE)

			null_modF <-suppressMessages(lmer(y~1+(1|perm_ID),dat[[i]], REML=FALSE))

			null_dat <- list(
				N = nrow(dat[[i]]),
				N_ID = length(unique(dat[[i]][,ID])),
				y = dat[[i]][,y],
				ID = perm_ID
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

		}))
		cat(i, " ")
		list(param = c(ICC=ICC, N_group=N_group, N_within=N_within), data=dat[[i]],actual=actual,null=null)

	},mc.cores=mc.cores)
	return(dist_M)
}



# ICC <- c(0,0.2,0.4,0.8)
# N_cluster <- c(20,40,80,160)
# N_within <- c(2,4,8,16)
system.time({
for(ICC in 0.2){
	for(N_group in c(20,40,80,160)){
		for(N_within in c(2,4,8,16)){
			cat("\n",ICC,N_group,N_within,"\n")
			sim_dat <- perm_sim(N_pop=100, N_perm=100, ICC=ICC, N_group=N_group, N_within=N_within)
			save(sim_dat, file=paste0(wd,"/Data/Intermediate/sims_",ICC,"_",N_group,"_",N_within,".Rdata"))
		}
	}
}
})

# system.time({dat1 <- perm_sim(N_pop=100, N_perm=100, ICC=0.2, N_group=20, N_within=2)})




c(freq=sum(sapply(dat1, function(x) 1-mean(x$actual[1]>x$null[,1]))>0.05),
mean=sum(sapply(dat1, function(x) 1-mean(x$actual[2]>x$null[,2]))>0.05),
mode=sum(1-sapply(dat1, function(x) mean(x$actual[3]>x$null[,3]))>0.05)
)

par(mfrow=c(3,1))
hist(dat1[[1]][["null"]][,1], xlim=c(0,1), breaks=30)
hist(dat1[[1]][["null"]][,2], xlim=c(0,1), breaks=30)
hist(dat1[[1]][["null"]][,3], xlim=c(0,1), breaks=30)


hist(dist_M[[3]][["null"]][,1], xlim=c(0,1), breaks=30)
hist(dist_M[[3]][["null"]][,2], xlim=c(0,1), breaks=30)

actual<-t(sapply(dat1, function(x) x$actual))
par(mfrow=c(3,1))
hist(actual[,1], xlim=c(0,1), breaks=50); abline(v=0.3, col="red",lwd=3)
hist(actual[,2], xlim=c(0,1), breaks=50); abline(v=0.3, col="red",lwd=3)
hist(actual[,3], xlim=c(0,1), breaks=50); abline(v=0.3, col="red",lwd=3)
colMeans(actual)



# prior <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=0.002)))
prior_PE <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alphaV=1000)))
# prior_flat <-list(R = list(V = 1, nu = 0), G=list(G1 = list(V = 1, nu = 0)))



dat<-dat[[1]]
m0 <- brm(y~1+(1|ID),dat, iter = 5000, warmup = 1000,save_pars = save_pars(all = TRUE), chains=1,backend="cmdstanr")

m1 <- brm(y~1,dat,chains = 1, iter = 5000, warmup = 1000,save_pars = save_pars(all = TRUE),backend="cmdstanr")

bayesfactor_models(m0, denominator = m1)

x <- seq(0,1,length.out=10000)
plot(x,dcauchy(x,0,2), ylim=c(0,max(dcauchy(x,0,2))))

{
par(mfcol=c(3,2))
hist(dist_M[1,], xlim=range(dist_M[1:3,]), breaks=30)
abline(v=0.3)
hist(dist_M[2,], xlim=range(dist_M[1:3,]), breaks=30)
abline(v=0.3)
hist(dist_M[3,], xlim=range(dist_M[1:3,]), breaks=30)
abline(v=0.3)
plot(dist_M[1,],dist_M[2,]); abline(0,1)
plot(dist_M[3,],dist_M[2,]); abline(0,1)
hist(dist_M[4,], breaks=30)
}
rowMeans(dist_M)
plot(dist_M[3,],dist_M[4,])




iw_prior <- function(v,V,nu)  MCMCpack::dinvgamma(v, shape = nu/2,scale =nu*V/2)
pe_prior <- function(v,V,nu,alpha.mu,alpha.V) df(v/alpha.V, df1 = 1, df2 = nu, ncp = (alpha.mu^2)/alpha.V)
pe_prior_sd <- function(v,V,nu,alpha.mu,alpha.V) 2 * dt(sqrt(v)/sqrt(alpha.V), df = nu, ncp = alpha.mu/sqrt(alpha.V))

x<-seq(min(dist_M),max(dist_M),length.out=10000)
plot(x,iw_prior(x,V=1,nu=0.002), type="l")
plot(x,pe_prior(x,V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), type="l")
plot(x,pe_prior_sd(x,V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), type="l")



mod<-lmer(y~1+(1|ID),dat)
var_out <- function(mod) c(as.numeric(summary(mod)$varcor),as.numeric(summary(mod)$sigma)^2)
var_out(mod)
var_outM <- function(mod) summary(mod)$Gcov[1]



nullmodelsVCV <- apply(1:100,function(i){
	dat$nully <- sample(dat$y, nrow(dat))
    nullmodel <- lmer(nully~1+(1|ID),dat)
     var_out(nullmodel)
})

hist(nullmodelsVCV[1,], breaks=100)
abline(v=var_out(mod)[1], col="red")
var_out(mod)[1]
1-sum(var_out(mod)[1]>nullmodelsVCV[1,])/100

mod<-lmer(y~1+(1|ID),dat)
mod1<-lm(y~1,dat)

bayesfactor_models(mod, denominator = mod1)

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

dist_M <- mclapply(1:100,function(i){
    mod <- MCMCglmm(y~1, random=~ID,data=dat, verbose=FALSE, nitt=13000*5, thin=10*5, burnin=3000*5)
    c(mean(mod$VCV[,1]), posterior.mode(mod$VCV[,1]))
    print(i)
})

par(mfrow=c(2,1))
hist(dist_M[1,], xlim=range(dist_M))
abline(v=0.3)
hist(dist_M[2,], xlim=range(dist_M))
abline(v=0.3)


# intercept only model
m0 <- brm(y~1+(1|ID),dat,chains = 10, iter = 5000, warmup = 1000,save_pars = save_pars(all = TRUE))
m1 <- brm(y~1,dat,chains = 10, iter = 5000, warmup = 1000,save_pars = save_pars(all = TRUE))

bayesfactor_models(m0, denominator = m1)




p_latent <- rnorm(10000,0,1)
p_exp <- pnorm(p_latent) 
hist(p_exp)

