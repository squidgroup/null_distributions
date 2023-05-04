rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

convert_func <- function(x){
	param <- if(length(x$param)==4) x$param[2:4] else x$param
	# param <-  x$param[2:4]
	names(param) <-c("ICC", "N_group", "N_within")
	c(param,
		median = as.numeric(x$actual["median"]),
		median_p = p_func(x$actual["median"],x$boot[,"median"]),
		median_bias = as.numeric(x$actual["median"] - param["ICC"]),
		median_rel_bias = as.numeric((x$actual["median"] - param["ICC"])/param["ICC"]),
		mean = as.numeric(x$actual["mean"]),
		mean_bias = as.numeric(x$actual["mean"] - param["ICC"]),
		mean_rel_bias = as.numeric((x$actual["mean"] - param["ICC"])/param["ICC"]),
		mode1 = as.numeric(x$actual["mode1"]),		
		mode1_bias = as.numeric(x$actual["mode1"] - param["ICC"]),
		mode1_rel_bias = as.numeric((x$actual["mode1"] - param["ICC"])/param["ICC"]),
		mode0.1 = as.numeric(x$actual["mode0.1"]),
		mode0.1_bias = as.numeric(x$actual["mode0.1"] - param["ICC"]),
		mode0.1_rel_bias = as.numeric((x$actual["mode0.1"] - param["ICC"])/param["ICC"])
	)
}

# files_p <- list.files(paste0(wd,"Data/Intermediate"))
# results_p <- files_p[grep("gaus_perm",files_p)]
# list_names_p <- gsub("gaus_perm_|.Rdata","",results_p)

# p_perm<-as.data.frame(do.call(rbind,lapply(results_p, function(i){
#   load(paste0(wd,"Data/Intermediate/",i))
#   t(sapply(perm_out, convert_func))
# })))

# head(p_perm)

files <- list.files(paste0(wd,"Data/Intermediate"))

results_gaus <- files[grep("gaus_boot",files)]
p_gaus<-as.data.frame(do.call(rbind,lapply(results_gaus, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  boot_out<-lapply(boot_out, function(x) {
  	x$param <- as.numeric(strsplit(gsub("gaus_boot_|.Rdata","",i),"_")[[1]])
		x
	})
  t(sapply(boot_out, convert_func))
})))

results_bern <- files[grep("bern_boot",files)]
p_bern<-as.data.frame(do.call(rbind,lapply(results_bern, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  # print(length(boot_out))
  t(sapply(boot_out[601:110], convert_func))
})))
#   load(paste0(wd,"Data/Intermediate/pois_boot_0.1.Rdata"))
# boot_out<-boot_out[201:400]
# save(boot_out,file=paste0(wd,"Data/Intermediate/pois_boot_0.1.Rdata"))


results_pois <- files[grep("pois_boot",files)]
p_pois<-as.data.frame(do.call(rbind,lapply(results_pois, function(i){
	## adjust ICC by total latent var
  load(paste0(wd,"Data/Intermediate/",i))
    boot_out<-lapply(boot_out, function(x) {
  	x$param["ICC"] <- x$param["ICC"]*0.2
		x
	})
  t(sapply(boot_out, convert_func))
})))


# i<-results_gaus[1]
# x<-boot_out[[1]]
# names(boot_out)
# convert_func(x)


# power_dat <- aggregate(median_p~N_group+ N_within+ICC,p_perm,function(x) mean(x<0.05))

# bias_dat <- aggregate(cbind(median_bias,median_rel_bias,mode_bias,mode_rel_bias)~N_group+ N_within+ICC,p_perm, mean)

# # plot(power_dat$median_p,bias_dat$median_bias)
# plot(power_dat$median_p,bias_dat$median_rel_bias)

# plot(power_dat$median_p,bias_dat$mode_rel_bias)
# # plot(power_dat$median_p,bias_dat$mode_bias)


## plot GLMM results on here - would help show its a power issue




# load(paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))
# load( file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))

# p_alt <- as.data.frame(t(sapply(results0.2, convert_func)))

# p_null <- as.data.frame(t(sapply(results0, convert_func
# 		)))
# head(p_null)



# load(paste0(wd,"Data/Intermediate/fay_results.Rdata"))
# results_sum <- t(sapply(results,colMeans))



p_bern$model <- "bern"
p_pois$model <- "pois"
p_gaus$model <- "gaus"
all <- rbind(p_bern,p_gaus,p_pois)

power_dat <- aggregate(median_p~N_group+ N_within+ICC+model,all,function(x) mean(x<0.05))
bias_dat <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, mean)
bias_CIs <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, function(x) 1.96*sd(x)/sqrt(length(x)))

prec_dat <- aggregate(cbind(median,mean,mode1,mode0.1)~N_group+ N_within+ICC+model,all,function(x) 1/sd(x))
prec_rel_dat <- aggregate(cbind(median,mean,mode1,mode0.1)~N_group+ N_within+ICC+model,all,function(x) mean(x)/sd(x))


plot_function <- function(metric){
	ylim <- range(bias_dat[,c("mean_rel_bias","median_rel_bias","mode1_rel_bias","mode0.1_rel_bias")], finite=TRUE)
	code <- paste0(metric,"_rel_bias")
	plot(power_dat$median_p,bias_dat[,code], pch=19, col=c(2,1,4)[as.factor(power_dat$model)], ylim=ylim, xlab="Power", ylab=paste(metric,"relative bias"))
	abline(h=0, col="grey")
	arrows(power_dat$median_p,bias_dat[,code] + bias_CIs[,code],power_dat$median_p,bias_dat[,code] - bias_CIs[,code], angle=90,code=3, length=0.05, col=c(2,1,4)[as.factor(power_dat$model)])	
}

plot_function2 <- function(metric){
	ylim <- range(bias_dat[,c("mean_rel_bias","median_rel_bias","mode1_rel_bias","mode0.1_rel_bias")], finite=TRUE)
	code <- paste0(metric,"_rel_bias")
	plot(prec_rel_dat[,metric],bias_dat[,code], pch=19, col=c(2,1,4)[as.factor(prec_rel_dat$model)], ylim=ylim, xlim=c(0,5), xlab="Relative Precision", ylab=paste(metric,"relative bias"))
	abline(h=0, col="grey")
	arrows(prec_rel_dat[,metric],bias_dat[,code] + bias_CIs[,code],prec_rel_dat[,metric],bias_dat[,code] - bias_CIs[,code], angle=90,code=3, length=0.05, col=c(2,1,4)[as.factor(prec_rel_dat$model)])	
}


{
par(mfrow=c(2,2), mar=c(5,5,1,1))

plot_function("mean")

plot_function("median")

plot_function("mode1")

plot_function("mode0.1")
legend("topright", c("Gaussian","Poisson","Bernoulli"),pch=19, col=c(1,4,2))
}

{
par(mfrow=c(2,2), mar=c(5,5,1,1))

plot_function2("mean")

plot_function2("median")

plot_function2("mode1")

plot_function2("mode0.1")
legend("topright", c("Gaussian","Poisson","Bernoulli"),pch=19, col=c(1,4,2))
}


	plot(power_dat$median_p,prec_dat[,"median"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")

	plot(power_dat$median_p,prec_rel_dat[,"median"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")

	plot(power_dat$median_p,prec_rel_dat[,"mean"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")

	plot(power_dat$median_p,prec_rel_dat[,"mode1"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")

	plot(power_dat$median_p,prec_rel_dat[,"mode0.1"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")


