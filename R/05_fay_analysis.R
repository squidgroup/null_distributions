
rm(list=ls())
library(squidSim)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

run_sims<- FALSE

if(run_sims){
	sim_dat <- simulate_population(
		data_structure = make_structure("individual(4000) + year(10)"),
		n_response=2,
		response_names=c("mortality","reproduction"),
		link="logit",
		family="binomial",
		parameters = list(
			intercept=c(qlogis(0.5),qlogis(0.7)),
			individual=list(
				vcorr=matrix(c(0.2^2,0.6,0.6,0.2^2),ncol=2)
			),
			residual=list(vcov=c(0,0))
		),
		n_pop=100,
		sample_type="nested",
		sample_param=as.matrix(expand.grid(individual=c(250,500,1000,2000,4000),observation=c(10))) 
	)
	# dat <- get_population_data(sim_dat, list=TRUE)

	surival_sampling <- function(dat, ID_name, mortality_name){
	  do.call(rbind,lapply(split(dat,dat[,ID_name]), function(x){
		  if(1 %in% x[,mortality_name]) {
		 		x[1:which(cumsum(x[,mortality_name])==1)[1],]
			}else {x}
		}))
	}

	results<-list()

	for (j in 1:nrow(sim_dat$sample_param)){
		dat <- get_sample_data(sim_dat, sample_set = j, list=TRUE)
		param <- sim_dat$sample_param[j,]
		cat("\n set",j,"\n")

		surv_dat<-lapply(dat,surival_sampling,ID_name="individual",mortality_name="mortality")

		GLMM_stan <- stan_model(file = paste0(wd,"stan/simple_GLMM.stan"))

		results[[j]] <- do.call(rbind,mclapply(1:length(surv_dat),function(i){
			surv_dat[[i]]$individual <- as.numeric(as.factor(surv_dat[[i]]$individual))

			stan_dat <- list(
				N = nrow(surv_dat[[i]]),
				N_ID = length(unique(surv_dat[[i]][,"individual"])),
				y = surv_dat[[i]][,"reproduction"],
				ID = surv_dat[[i]][,"individual"])

			stan_mod <- sampling(GLMM_stan, data=stan_dat, chains=1,iter=7500, warmup=2500, pars=c("sigma2_ID"), refresh=0)

			cat(i, " ")

			c(n=as.numeric(param[1]),stan_out(stan_mod))
				

		},mc.cores=8))
	}

	save(results, file=paste0(wd,"Data/Intermediate/fay_results.Rdata"))
}
