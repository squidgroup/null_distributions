
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))
pops <- 101:300
sets <- expand.grid(ICC=c(0,0.1,0.2,0.4),N_group=c(20,40,80),N_within=c(2,4))
N_boot=100

for(j in 4:nrow(sets)){
	cat("\n set",j,":",as.matrix(sets[j,]), "\n")

	file_end <- paste0(sets[j,"ICC"],"_",sets[j,"N_group"],"_",sets[j,"N_within"],".Rdata")
	
	load(paste0(wd,"Data/Intermediate/gaus_sims_",file_end))

	boot_out <- mclapply(pops,function(i){
		set.seed(paste0(i,paste(sets[j,c(2,3,1)], collapse="")))
		cat(i, " ")

		bs_dat_all<-simulate_population(
			data_structure= make_structure(paste0("ID(",out[[i]]$param[3],")"),repeat_obs=out[[i]]$param[4]),
			parameters= list(residual=list(vcov=median(rowSums(out[[i]]$post)))),
			n_pop=N_boot
		)

		bootstrap <- do.call(rbind,lapply(get_population_data(bs_dat_all, list=TRUE),function(j){
			
			c(gaussian_mods(j)$summary)

		}))

		list(
			param = out[[i]]$param,
			actual=out[[i]]$summary,
			bootstrap=bootstrap)

	},mc.cores=8)

  # is it exists already, join them together and save the two	
	file <- paste0("gaus_boot_",file_end)
	if(file %in% list.files(paste0(wd,"Data/Intermediate/"))){
		boot_out2 <- boot_out
		load(paste0(wd,"Data/Intermediate/",file))
		boot_out <- c(boot_out, boot_out2)
	}

	save(boot_out, file=paste0(wd,"Data/Intermediate/",file))
}




#	save(PB_out, file=paste0(wd,"Data/Intermediate/gaus_PB_",file_end))

