rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

convert_func <- function(x){
	param <- if(length(x$param)==4) x$param[2:4] else x$param
	names(param) <-c("ICC", "N_group", "N_within")
	c(param,
		median_p = p_func(x$actual["median"],x$perm[,"median"]),
		median_bias = as.numeric(x$actual["median"] - param["ICC"]),
		median_rel_bias = as.numeric((x$actual["median"] - param["ICC"])/param["ICC"]),
		mean_bias = as.numeric(x$actual["mean"] - param["ICC"]),
		mean_rel_bias = as.numeric((x$actual["mean"] - param["ICC"])/param["ICC"]),
		mode1_bias = as.numeric(x$actual["mode1"] - param["ICC"]),
		mode1_rel_bias = as.numeric((x$actual["mode1"] - param["ICC"])/param["ICC"]),
		mode0.1_bias = as.numeric(x$actual["mode0.1"] - param["ICC"]),
		mode0.1_rel_bias = as.numeric((x$actual["mode0.1"] - param["ICC"])/param["ICC"])
	)
}

files_p <- list.files(paste0(wd,"Data/Intermediate"))
results_p <- files_p[grep("gaus_perm",files_p)]
list_names_p <- gsub("gaus_perm_|.Rdata","",results_p)

p_perm<-as.data.frame(do.call(rbind,lapply(results_p, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(perm_out, convert_func))
})))


head(p_perm)

# power_dat <- aggregate(median_p~N_group+ N_within+ICC,p_perm,function(x) mean(x<0.05))

# bias_dat <- aggregate(cbind(median_bias,median_rel_bias,mode_bias,mode_rel_bias)~N_group+ N_within+ICC,p_perm, mean)

# # plot(power_dat$median_p,bias_dat$median_bias)
# plot(power_dat$median_p,bias_dat$median_rel_bias)

# plot(power_dat$median_p,bias_dat$mode_rel_bias)
# # plot(power_dat$median_p,bias_dat$mode_bias)


## plot GLMM results on here - would help show its a power issue




# load(paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))
load( file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))

p_alt <- as.data.frame(t(sapply(results0.2, convert_func)))

p_null <- as.data.frame(t(sapply(results0, convert_func
		)))
head(p_null)



load(paste0(wd,"Data/Intermediate/fay_results.Rdata"))
results_sum <- t(sapply(results,colMeans))



p_null$model <- p_alt$model <- "glmm"
p_perm$model <- "lmm"
all <- rbind(p_perm,p_alt,p_null)

power_dat <- aggregate(median_p~N_group+ N_within+ICC+model,all,function(x) mean(x<0.05))
bias_dat <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, mean)
bias_CIs <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, function(x) 1.96*sd(x)/sqrt(length(x)))


plot_function <- function(metric){
	ylim <- range(bias_dat[,c("mean_rel_bias","median_rel_bias","mode1_rel_bias","mode0.1_rel_bias")], finite=TRUE)
	code <- paste0(metric,"_rel_bias")
	plot(power_dat$median_p,bias_dat[,code], pch=19, col=c(2,1)[as.factor(power_dat$model)], ylim=ylim, xlab="Power")
	abline(h=0, col="grey")
	arrows(power_dat$median_p,bias_dat[,code] + bias_CIs[,code],power_dat$median_p,bias_dat[,code] - bias_CIs[,code], angle=90,code=3, length=0.05)
	
}


{
par(mfrow=c(2,2), mar=c(5,5,1,1))

plot_function("mean")

plot_function("median")

plot_function("mode1")

plot_function("mode0.1")
}

