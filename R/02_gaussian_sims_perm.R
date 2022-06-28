
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
pops <- 301:500
sets <- expand.grid(ICC=c(0,0.1,0.2,0.4),N_group=c(20,40,80),N_within=c(2,4))
N_perm=100

for(j in 1:nrow(sets)){
	cat("\n set",j,":",as.matrix(sets[j,]), "\n")

	file_end <- paste0(sets[j,"ICC"],"_",sets[j,"N_group"],"_",sets[j,"N_within"],".Rdata")
	
	load(paste0(wd,"Data/Intermediate/gaus_sims_",file_end))

	perm_out <- mclapply(pops,function(i){
		set.seed(paste0(i,paste(sets[j,c(2,3,1)], collapse="")))
		cat(i, " ")
		data <- out[[i]]$data

		perm <- do.call(rbind,lapply(1:N_perm,function(j){
			data[,"ID"] <- sample(data[,"ID"], replace=FALSE)
	  	c(gaussian_mods(data)$summary)
		}))

		list(
			param = out[[i]]$param,
			actual = out[[i]]$summary,
			perm = perm
		)

	},mc.cores=8)

  # is it exists already, join them together and save the two	
	file <- paste0("gaus_perm_",file_end)
	if(file %in% list.files(paste0(wd,"Data/Intermediate/"))){
		perm_out2 <- perm_out
		load(paste0(wd,"Data/Intermediate/",file))
		perm_out <- c(perm_out, perm_out2)
	}

	save(perm_out, file=paste0(wd,"Data/Intermediate/",file))
}




#	save(PB_out, file=paste0(wd,"Data/Intermediate/gaus_PB_",file_end))

