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
  t(sapply(boot_out, convert_func))
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


p_bern$model <- "bern"
p_pois$model <- "pois"
p_gaus$model <- "gaus"
all <- rbind(p_bern,p_gaus,p_pois)






power_dat <- aggregate(median_p~N_group+ N_within+ICC+model,all,function(x) mean(x<0.05))
bias_dat <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, mean)
bias_CIs <- aggregate(cbind(median_bias,median_rel_bias,mean_bias,mean_rel_bias,mode1_bias,mode1_rel_bias,mode0.1_bias,mode0.1_rel_bias)~N_group+ N_within+ICC+model,all, function(x) 1.96*sd(x)/sqrt(length(x)))


plot_function <- function(metric){
	# ylim <- range(bias_dat[,c("mean_rel_bias","median_rel_bias","mode1_rel_bias","mode0.1_rel_bias")], finite=TRUE)
	ylim <- c(-1,1.5)
	code <- paste0(metric,"_rel_bias")
	plot(power_dat$median_p,bias_dat[,code], pch=19, col=c(2,1,4)[as.factor(power_dat$model)], ylim=ylim, xlab="Power", ylab=paste(metric,"relative bias"))
	abline(h=0, col="grey")
	arrows(power_dat$median_p,bias_dat[,code] + bias_CIs[,code],power_dat$median_p,bias_dat[,code] - bias_CIs[,code], angle=90,code=3, length=0.05, col=c(2,1,4)[as.factor(power_dat$model)])	
}



setEPS()
pdf(paste0(wd,"Figures/Fig6_power_bias.pdf"), height=8, width=10)
{
par(mfrow=c(2,2), mar=c(4,4,1,1), cex.axis=1,cex.lab=1.25, mgp=c(2,0.5,0))

plot_function("mean")

plot_function("median")

plot_function("mode1")

plot_function("mode0.1")
legend("topright", c("Gaussian","Poisson","Bernoulli"),pch=19, col=c(1,4,2), cex=1.25)
}
dev.off()



## Relative Precision and relative bias

prec_dat <- aggregate(cbind(median,mean,mode1,mode0.1)~N_group+ N_within+ICC+model,all,function(x) 1/sd(x))
prec_rel_dat <- aggregate(cbind(median,mean,mode1,mode0.1)~N_group+ N_within+ICC+model,all,function(x) mean(x)/sd(x))


plot_function2 <- function(metric,...){
	# ylim <- range(bias_dat[,c("mean_rel_bias","median_rel_bias","mode1_rel_bias","mode0.1_rel_bias")], finite=TRUE)
	ylim <- c(-1,1.5)
	code <- paste0(metric,"_rel_bias")
	plot(prec_rel_dat[,metric],bias_dat[,code], pch=19, col=c(2,1,4)[as.factor(prec_rel_dat$model)], ylim=ylim, xlim=c(0,5), xlab="Relative Precision", ylab="Relative bias",...)
	abline(h=0, col="grey")
	arrows(prec_rel_dat[,metric],bias_dat[,code] + bias_CIs[,code],prec_rel_dat[,metric],bias_dat[,code] - bias_CIs[,code], angle=90,code=3, length=0.05, col=c(2,1,4)[as.factor(prec_rel_dat$model)])	
}


setEPS()
pdf(paste0(wd,"Figures/FigSM_precision_bias.pdf"), height=8, width=10)
{
par(mfrow=c(2,2), mar=c(3,3,3,1), cex.axis=1,cex.lab=1.25, mgp=c(2,0.5,0))
plot_function2("mean", main="Mean")

plot_function2("median", main="Median")

plot_function2("mode1", main="Mode-1")

plot_function2("mode0.1", main="Mode-0.1")
legend("topright", c("Gaussian","Poisson","Bernoulli"),pch=19, col=c(1,4,2))
}
dev.off()

	# plot(power_dat$median_p,prec_dat[,"median"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Power")

setEPS()
pdf(paste0(wd,"Figures/FigSM_precision_power.pdf"), height=8, width=10)
{
par(mfrow=c(2,2), mar=c(3,3,3,1), cex.axis=1,cex.lab=1.25, mgp=c(2,0.5,0))
	
	plot(power_dat$median_p,prec_rel_dat[,"mean"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Relative Precision",  ylab="Power", ylim=c(0,5), main="Mean")

	plot(power_dat$median_p,prec_rel_dat[,"median"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Relative Precision",  ylab="Power", ylim=c(0,5), main="Median")

	plot(power_dat$median_p,prec_rel_dat[,"mode1"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Relative Precision",  ylab="Power", ylim=c(0,5), main="Mode-1")

	plot(power_dat$median_p,prec_rel_dat[,"mode0.1"], pch=19, col=c(2,1,4)[as.factor(power_dat$model)],  xlab="Relative Precision",  ylab="Power", ylim=c(0,5), main="Mode-0.1")
}
dev.off()




##### p distribution


p_mean_dat <- aggregate(median_p~N_group+ N_within+ICC+model,all,mean)

p_var_dat <- aggregate(median_p~N_group+ N_within+ICC+model,all,var)


setEPS()
pdf(paste0(wd,"Figures/FigSM_p_dist.pdf"), height=5, width=10)

{
par(mfrow=c(1,2), mar=c(5,5,1,1), cex.axis=0.75, mgp=c(2,0.5,0))
plot(power_dat$median_p,p_mean_dat$median_p, pch=c(17,rep(19,7))[as.factor(power_dat$ICC)], col=c(2,1,4)[as.factor(power_dat$model)], ylab="Mean p-value", xlab="Power")

plot(power_dat$median_p,p_var_dat$median_p, pch=c(17,rep(19,7))[as.factor(power_dat$ICC)], col=c(2,1,4)[as.factor(power_dat$model)], ylab="Variance in p-values", xlab="Power")

legend("topright", c("Gaussian","Poisson","Bernoulli","ICC > 0","ICC = 0"),pch=c(rep(19,4),17), col=c(1,4,2,"grey","grey"))
}
dev.off()




# alpha <- 1
# mu = seq(0,alpha,0.01)
# power = pbeta(0.05,1,(alpha/mu -alpha))#1-(1-0.05)^(1/mu -1)
# lines(mu~power)



# plot(power_dat2$median_p,power_dat3$median_p, pch=19, col=c(2,1,4)[as.factor(power_dat$model)], ylab="Variance in p-values", xlab="Mean p-value")

# mu = seq(0,0.5,0.01)
#  V = (mu^2-mu^3)/(1+mu)
#  lines(V~mu)

# alpha <- 1
# V = mu
