
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))



perm_boot <- function(output, N_perm=100, N_boot=100){
#output<-out[[1]]
	data <- output$data
  ICC <- output$ICC
  N_group <- output$N_group
  N_within <- output$N_within

	perm <- do.call(rbind,lapply(1:N_perm,function(j){
		data[,"ID"] <- sample(data[,ID], replace=FALSE)

		c(gaussian_mods(data)$summary)
	}))

	bs_dat_all<-simulate_population(
		data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
		parameters= list(residual=list(vcov=median(rowSums(output$post)))),
		n_pop=N_boot
	)

	bootstrap <- do.call(rbind,lapply(get_population_data(bs_dat_all, list=TRUE),function(j){
		
		c(gaussian_mods(j)$summary)

	}))

	list(param = c(ICC=ICC, N_group=N_group, N_within=N_within),actual=output$summary,perm=perm, bootstrap=bootstrap)

}


LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))
pops <- 101:102
sets <- expand.grid(ICC=c(0,0.2,0.4),N_group=c(20,40,80,160),N_within=c(2,4))

for(j in 1:24){
	cat("\n set:",as.matrix(sets[j,]), "\n")

	file_end <- paste0(sets[j,"ICC"],"_",sets[j,"N_group"],"_",sets[j,"N_within"],".Rdata")
	
	load(paste0(wd,"Data/Intermediate/gaus_sims_",file_end))

	set.seed(paste0("2206",j,min(pops)))

	PB_out <- mclapply(pops,function(i){
		cat(i, " ")
		perm_boot(out[[i]])
	},mc.cores=8)

# is it exists already, join them together and save the two	
	file <- paste0("gaus_PB_",file_end)
	if(file %in% list.files(paste0(wd,"Data/Intermediate/"))){
		PB_out2 <- PB_out
		load(paste0(wd,"Data/Intermediate/",file))
		PB_out <- c(PB_out, PB_out2)
	}

	save(PB_out, file=paste0(wd,"Data/Intermediate/",file))
}




#	save(PB_out, file=paste0(wd,"Data/Intermediate/gaus_PB_",file_end))

