rm(list=ls())

library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("sims",files)]
list_names <- gsub("sims_|.Rdata","",results)

## get all results into a list
all <- list()
for(i in 1:length(results)){
	load(paste0(wd,"Data/Intermediate/",results[i]))
	all[[list_names[i]]] <- sim_dat
	rm("sim_dat")
}

power <- as.data.frame(t(sapply(all,function(z){
	c(
		z[[1]]$param,LRT=mean(sapply(z, function(x) x$actual["LRT"])< 0.05),
		freq=mean(sapply(z, function(x) p_func(x$actual["freq"],x$null[,"freq"]) < 0.05)),
		mean=mean(sapply(z, function(x) p_func(x$actual["mean"],x$null[,"mean"]) < 0.05)),
		mode=mean(sapply(z, function(x) p_func(x$actual["mode"],x$null[,"mode"]) < 0.05))
	)
})))
power <- power[order(power$ICC,power$N_group,power$N_within),]
power<-subset(power,ICC!=0.8 & N_within!=16)


par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(freq~mean,power, xlab="Bayesian Power", ylab="ML Power"); abline(0,1, col="grey")
points(LRT~mean,power, pch=19, col="red")
