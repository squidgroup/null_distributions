
rm(list=ls())

devtools::load_all("~/github/squid/R")

library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))


set.seed(1250)
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

		actual<-stan_out(stan_mod)

	null <- do.call(rbind,mclapply(1:1000,function(j){
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


squid_dat2<-simulate_population(
		data_structure= make_structure(paste0("ID(80)"),repeat_obs=2),
		parameters= list( ID=list(vcov=0.4), residual=list(vcov=0.6)),
		N_pop=1
		)
dat2 <- get_population_data(squid_dat2)

stan_dat2 <- list(
			N = nrow(dat2),
			N_ID = length(unique(dat2[,ID])),
			y = dat2[,y],
			ID = dat2[,ID])

stan_mod2 <- sampling(LMM_stan, data=stan_dat2, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)

		actual2<-stan_out(stan_mod2)

	null2 <- do.call(rbind,mclapply(1:1000,function(j){
			perm_ID <- sample(dat2[,ID], replace=FALSE)

			null_modF <-suppressMessages(lmer(y~1+(1|perm_ID),dat2, REML=FALSE))

			null_dat <- list(
				N = nrow(dat2),
				N_ID = length(unique(dat2[,ID])),
				y = dat2[,y],
				ID = perm_ID
			)

			null_mod <- sampling(LMM_stan, data=null_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
			c(freq=as.numeric(summary(null_modF)$varcor),stan_out(null_mod))

		}))

save(actual,actual2,null,null2,file=paste0(wd,"Data/Intermediate/figure1_data.Rdata"))


p_perm <- p_func(actual["mean"],null[,"mean"])

setEPS()
pdf(paste0(wd,"Figures/Fig1.pdf"), height=6, width=6)

{
par(mfcol=c(2,2), mar=c(5,1,1,1))
hist(extract(stan_mod)$sigma2_ID, main="",xlim=c(0,1),freq=FALSE, breaks=seq(0,2,length.out=101),ylim=c(0,6), xlab="Posterior Samples", yaxt="n", ylab="")
abline(v=actual[c("mean","mode")], col=c("red","blue"),lwd=2)
text(0.8,4,print_func(actual))
mtext("a)",2,padj=-7, las=1)


hist(null[,"mean"], main="",xlim=c(0,1),freq=FALSE, breaks=seq(0,2,length.out=101),ylim=c(0,16), xlab="Permuted Posterior Means", yaxt="n", ylab="")
abline(v=actual["mean"], col="red",lwd=2)
text(0.8,10,bquote(P[perm]==.(p_perm)))
mtext("c)",2,padj=-7, las=1)


hist(extract(stan_mod2)$sigma2_ID, main="",xlim=c(0,1),freq=FALSE, breaks=seq(0,2,length.out=101),ylim=c(0,6), xlab="Posterior Samples", yaxt="n", ylab="")
abline(v=actual2[c("mean","mode")], col=c("red","blue"),lwd=2)

text(0.8,4,print_func(actual2))
mtext("b)",2,padj=-7, las=1)


hist(null2[,"mean"], main="",xlim=c(0,1),freq=FALSE, breaks=seq(0,2,length.out=101),ylim=c(0,16), xlab="Permuted Posterior Means", yaxt="n", ylab="")
abline(v=actual2["mean"], col="red",lwd=2)
p_perm2 <- p_func(actual2["mean"],null2[,"mean"])
 # text(0.8,10,bquote(P[perm][.(p_perm)]))
  text(0.8,10,bquote(P[perm]==.(p_perm2)))
mtext("d)",2,padj=-7, las=1)

}
dev.off()






# # prior <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=0.002)))
# prior_PE <-list(R=list(V=diag(1), nu=0.002),G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alphaV=1000)))
# # prior_flat <-list(R = list(V = 1, nu = 0), G=list(G1 = list(V = 1, nu = 0)))
# modM<-MCMCglmm(y~1, random=~ID,data=dat)
# # modM1<-MCMCglmm(y~1,data=dat)
# summary(modM)$Gcov[1]
# plot(modM)

# nullmodelsVCV_M <- sapply(1:100,function(i){
# 	dat2<-dat
# 	dat2$nully <- sample(dat2$y, nrow(dat2))
#     nullmodel <- MCMCglmm(nully~1, random=~ID,data=dat2)
#      var_outM(nullmodel)
# })

# nullmodelsVCV_M
# hist(nullmodelsVCV_M, breaks=100)
# abline(v=summary(modM)$Gcov[1], col="red")

