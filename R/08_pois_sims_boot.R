
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

pois_stan <- stan_model(file = paste0(wd,"stan/simple_pois_GLMM.stan"))
pops <- 1:500
ICCs <- c(0)
N_boot=100


# j<-0
# i=1
for(j in ICCs){
	cat("\n ICC:",j,"\n")

	load(paste0(wd,"Data/Intermediate/pois_sims_",j,".Rdata"))

	boot_out <- mclapply(pops,function(i){
		set.seed(paste0(i,j, collapse=""))
		cat(i, " ")

		bs_dat_all<-simulate_population(
			data_structure= make_structure("ID(80)",repeat_obs=2),
			parameters= list( 
				intercept= median(out[[i]]$post[,"beta_0"]),			
				residual=list(vcov=median(rowSums(out[[i]]$post[,c("sigma2_ID","sigma2_E")])))
			),
			family="poisson", link="log",
			n_pop=N_boot
		)

		bootstrap <- do.call(rbind,lapply(get_population_data(bs_dat_all, list=TRUE),function(j){	
			c(pois_mods(j)$summary)
		}))

		list(
			param = out[[i]]$param,
			actual=out[[i]]$summary,
			bootstrap=bootstrap)

	},mc.cores=8)

  # is it exists already, join them together and save the two	
	file <- paste0("pois_boot_",j,".Rdata")
	if(file %in% list.files(paste0(wd,"Data/Intermediate/"))){
		boot_out2 <- boot_out
		load(paste0(wd,"Data/Intermediate/",file))
		boot_out <- c(boot_out, boot_out2)
	}

	save(boot_out, file=paste0(wd,"Data/Intermediate/",file))
	cat("\n")
}




#	save(PB_out, file=paste0(wd,"Data/Intermediate/gaus_PB_",file_end))

